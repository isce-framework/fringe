#include "sbas.hpp"
#include <armadillo>
#include <algorithm>
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_port.h"
#include "fringe/fringe_common.hpp"

void process(sbasOptions *opts)
{
    int blockSize = 128;

    //Parse the pairs metadata to extract dates
    //And to build design matrix
    std::vector<std::string> dates;
    std::vector<Scene> scenes;
    Inverter sbasInv;

    //Connectivity matrix block
    {
        std::vector<std::vector<std::string>> pairs;
        for(int ii=0; ii< opts->nPairs; ii++)
        {
            std::vector<std::string> pair;
            pair.push_back(opts->pairs[ii].masterDate);
            pair.push_back(opts->pairs[ii].slaveDate);
            pairs.push_back(pair);

            //Check master date
            if (std::find(opts->dates.begin(), opts->dates.end(), pair[0]) == opts->dates.end())
            {
                opts->dates.push_back( pair[0]);
            }

            //Check slave date
            if (std::find(opts->dates.begin(), opts->dates.end(), pair[1]) == opts->dates.end())
            {
                opts->dates.push_back( pair[1]);
            }
        }

        //Sort the date
        std::sort(opts->dates.begin(), opts->dates.end());
        opts->nSAR = opts->dates.size();

        //Setup matrix sizes
        //Assuming quadratic for now
        sbasInv.setup(opts->nPairs, opts->nSAR,2);


        //Convert to dates relative to ref/first date
        //And convert time to years
        for(int ii=0; ii< opts->dates.size(); ii++)
        {
            sbasInv.epochs[ii] = daysSince1900(opts->dates[ii]);
        }

        if (opts->refDate.empty())
        {
            sbasInv.refIndex = 0;
            opts->refDate = opts->dates[0];
        }
        else
        {
            ptrdiff_t posr = std::find(opts->dates.begin(), opts->dates.end(), opts->refDate) - opts->dates.begin();
            if (posr >= opts->dates.size())
                throw std::invalid_argument("Invalid reference date");

            sbasInv.refIndex = posr;
        }
        sbasInv.epochs -= sbasInv.epochs[sbasInv.refIndex];
        sbasInv.epochs /= 365.25;

        //Populate incidence matrix
        for(int ii=0; ii<opts->nPairs; ii++)
        {
            ptrdiff_t posm = std::find(opts->dates.begin(), opts->dates.end(), pairs[ii][0]) - opts->dates.begin();
            ptrdiff_t poss = std::find(opts->dates.begin(), opts->dates.end(), pairs[ii][1]) - opts->dates.begin();
            sbasInv.Jmat(ii, posm) = 1.0;
            sbasInv.Jmat(ii, poss) = -1.0;
            opts->pairs[ii].deltaT = sbasInv.epochs[posm] - sbasInv.epochs[poss];
        }

        std::cout << "Number of unique dates = " << opts->nSAR << std::endl;

    } //End of connectivity matrix block

    //Ensure that we have Bperps for all the dates identified.
    if (opts->estimateDEMError)
    {
        scenes.resize(opts->nSAR);
        std::cout << "Input Scene list size: " << opts->scenes.size() << "\n";

        for(int ii=0; ii<opts->nSAR; ii++)
        {
            int ind=-1;
            for(int jj=0; jj < opts->scenes.size(); jj++)
            {
                if (opts->dates[ii].compare(opts->scenes[jj].date) == 0)
                {
                    ind = jj;
                    break;
                }
            }
            if (ind < 0)
            {
                std::cout << "Could not find entry for date " << opts->dates[ii] << "\n";
                throw std::invalid_argument("Could not find bperp for given date");
            }
            scenes[ii] = opts->scenes[ind];
        }

        std::cout << "Scene entries for all dates found \n";
    }
    


    //Get number of threads available
    int nThreads = numberOfThreads();
    std::cout << "Number of threads = " << nThreads << "\n";

    std::cout << "Ref date: " << opts->refDate << "\n";
    std::cout << "Limits: " << opts->dates[0] << " "
              << opts->dates[ opts->nSAR-1] << "\n";

    std::cout << "Epoch limits (yrs): " << sbasInv.epochs.min() << "  " << sbasInv.epochs.max() << "\n"; 
    

    //Prepare connectivity matrix and check rank
    sbasInv.prepare();

    
    //Negative checks for bboxes
    if (opts->refbox.size() != 4 )
    {
        std::cout << "Reference box should be provided as 4 elements \n"
                  << "[xoff, yoff, xsize, ysize] \n";

        throw std::invalid_argument("Refbox is not size 4");
    }   
    
    if (std::any_of(opts->refbox.begin(), opts->refbox.end(), [](int i){return i<0;}))
    {
        throw std::invalid_argument("Refbox includes negative values");
    }

    if (opts->bbox.size() != 4 )
    {
        std::cout << "Bounding reference box should be provided as 4 elements \n"
                  << "[xoff, yoff, xsize, ysize] \n";

        throw std::invalid_argument("Bbox is not size 4");
    }   
    
    if (std::any_of(opts->bbox.begin(), opts->bbox.end(), [](int i){return i<0;}))
    {
        throw std::invalid_argument("Bbox includes negative values");
    }


    //Start GDAL
    GDALAllRegister();

    //File reading and handling
    std::vector<GDALDataset*> unwDataset(opts->nPairs);
    std::vector<GDALDataset*> bperpDataset(opts->nSAR);
    GDALDataset* incDataset;
    GDALDataset* cohDataset;

    int inCols, inRows;
    
    std::cout << "Reading files and setting up references \n";
    for(int ii=0; ii < opts->nPairs; ii++)
    {
        unwDataset[ii] = reinterpret_cast<GDALDataset*> (GDALOpen( opts->pairs[ii].ifgName.c_str(), GA_ReadOnly));
        if (unwDataset[ii] == NULL)
        {
            std::cout << "Cannot open file { " << opts->pairs[ii].ifgName << "} for reading. \n";
            for (int jj=0; jj < ii; jj++)
            {
                GDALClose( unwDataset[jj]);
            }
            GDALDestroyDriverManager();

            throw std::invalid_argument("Invalid unwrapped file name encountered");
        }

        //Gather dimensions of the first image
        //Assumption is that all images are the same size
        if (ii == 0)
        {
            inCols = unwDataset[0]->GetRasterXSize();
            inRows = unwDataset[0]->GetRasterYSize();

            std::cout << "Input Raster Size: " << inCols << "P x " << inRows << "L\n";

            //Compare against refbox
            if ((opts->refbox[2] != 0) && (opts->refbox[3] != 0))
            {
                if ((opts->refbox[0] + opts->refbox[2]) >= inCols)
                    throw std::invalid_argument("Refbox exceeds limits in X");

                if ((opts->refbox[1] + opts->refbox[3]) >= inRows)
                    throw std::invalid_argument("Refbox exceeds limits in Y");
            }


            //Compare against bbox
            if (std::any_of(opts->bbox.begin(), opts->bbox.end(), [](int i){return i!=0;}))
            {
                
                if ((opts->bbox[0] ==0) && (opts->bbox[2] == 0))
                {
                    opts->bbox[2] = inCols;
                }
                else
                {
                    if ((opts->bbox[0] + opts->bbox[2]) >= inCols)
                        throw std::invalid_argument("Bbox exceeds limits in X");

                    if (opts->bbox[2] == 0)
                        opts->bbox[2] = inCols - opts->bbox[0];
                }


                if ((opts->bbox[1] == 0) && (opts->bbox[3] == 0))
                {
                    opts->bbox[3] = inRows;
                }
                else
                {
                    if ((opts->bbox[1] + opts->bbox[3]) >= inRows)
                        throw std::invalid_argument("Bbox exceeds limits in Y");

                    if (opts->bbox[3] == 0)
                        opts->bbox[3] = inRows - opts->bbox[1];
                }
            }
            else
            {
                //Processing full image
                opts->bbox = {0,0,inCols,inRows};
            }

        }

        if (!opts->pairs[ii].cohName.empty())
        {
            cohDataset = reinterpret_cast<GDALDataset*> (GDALOpen( opts->pairs[ii].cohName.c_str(), GA_ReadOnly));

            if (cohDataset == NULL)
            {
                std::cout << "Cannot open file { " << opts->pairs[ii].cohName << " } for reading. \n";
                for(int jj=0; jj <= ii; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }
                GDALDestroyDriverManager();

                throw std::invalid_argument("Invalid coherence file name encountered");
            }
        }
        else
        {
            cohDataset = NULL;
        }

        
        //Check that input datasets are same size
        {
            if ((unwDataset[ii]->GetRasterXSize() != inCols) || 
                (unwDataset[ii]->GetRasterYSize() != inRows))
            {
                opts->pairs[ii].print();

                for(int jj=0; jj<=ii; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }
                GDALClose(cohDataset);
                GDALDestroyDriverManager();
                throw std::invalid_argument("Unwrap file size mismatch");
            }

            if (cohDataset != NULL)
            {
                if ((cohDataset->GetRasterXSize() != inCols) || 
                    (cohDataset->GetRasterYSize() != inRows))
                {
                    opts->pairs[ii].print();

                    for(int jj=0; jj<=ii; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }
                    GDALClose(cohDataset);
                    GDALDestroyDriverManager();
                    throw std::invalid_argument("Coherence file size mismatch");
                }
            }
        }


        //Referencing of the datasets
        {
            if ((opts->refbox[2] == 0) || (opts->refbox[3] == 0))
            {
                std::cout << "Either one or both of xsize/ysize set to 0 for refbox \n";
                std::cout << "Assuming the data is pre-referenced \n";

                opts->pairs[ii].referenceOffset = 0.0;
            }
            else
            {

                arma::fmat unwref(opts->refbox[2], opts->refbox[3]);
                arma::fmat cohref(opts->refbox[2], opts->refbox[3]);

                if (ii==0) //Report only once
                {
                    std::cout << "Reference region location: " << "\n";
                    std::cout << "Input Lines: " << opts->refbox[1] << " - " << opts->refbox[1] + opts->refbox[3] << "\n";
                    std::cout << "Input Pixels: " << opts->refbox[0] << " - " << opts->refbox[0] + opts->refbox[2] << "\n";
                    std::cout << "Output Lines: " << opts->refbox[1] - opts->bbox[1] << " - " << opts->refbox[1] + opts->refbox[3] - opts->bbox[1] << "\n";
                    std::cout << "Output Pixels: " << opts->refbox[0] - opts->bbox[0] << " - " << opts->refbox[0] + opts->refbox[2] - opts->bbox[0] << "\n";
                }

                int status = unwDataset[ii]->GetRasterBand( unwDataset[ii]->GetRasterCount())->RasterIO( GF_Read, opts->refbox[0], opts->refbox[1],
                    opts->refbox[2], opts->refbox[3],
                    (void*) unwref.memptr(),
                    opts->refbox[2], opts->refbox[3], GDT_Float32,
                    sizeof(float), sizeof(float)*opts->refbox[2], NULL);

                if (status != 0)
                {
                    std::cout << "Error reading ref from image { " << opts->pairs[ii].ifgName << "\n";

                    for(int jj=0; jj<= ii; jj++)
                    {
                        GDALClose(unwDataset[ii]);
                    }
                    if (cohDataset!=NULL) GDALClose(cohDataset);

                    throw std::invalid_argument("Error reading reference box observations");
                }

                if (cohDataset != NULL)
                {

                    status = cohDataset->GetRasterBand( cohDataset->GetRasterCount())->RasterIO( GF_Read, opts->refbox[0], opts->refbox[1],
                      opts->refbox[2], opts->refbox[3],
                      cohref.memptr(),
                      opts->refbox[2], opts->refbox[3], GDT_Float32,
                      sizeof(float), sizeof(float)*opts->refbox[2], NULL);

                    
                    if (status != 0)
                    {
                        std::cout << "Error reading ref from image { " << opts->pairs[ii].cohName << " } \n";
                        for(int jj=0; jj <= ii; jj++)
                        {
                            GDALClose(unwDataset[jj]);
                        }
                        GDALClose(cohDataset);

                        throw std::invalid_argument("Error reading reference box coherence");
                    }

                    //Close the coherence dataset here since we are not going to use it anymore
                    GDALClose(cohDataset);
                    cohDataset = NULL;
                 }
                else
                {
                    //Make it slightly greater than threshold
                    cohref.fill(opts->pairs[ii].threshold + 0.01f);
                }


                //Estimate coherence mean
                arma::uvec validpix = arma::find( cohref > opts->pairs[ii].threshold);
                
                if (validpix.n_elem == 0)
                {
                    std::cout << "Referencing of pair failed \n";
                    std::cout << "Ifgname : " << opts->pairs[ii].ifgName << "\n";
                    std::cout << "Cohname : " << opts->pairs[ii].cohName << "\n";
                    throw std::invalid_argument("Referencing of pair failed");
                }

                opts->pairs[ii].referenceOffset = arma::mean(unwref(validpix));
            }
        }

        GDALTermProgress( ii/ (1.0 * opts->nPairs), NULL, NULL);

    }       

    std::cout << "\nCompleted setting up references. \n"; 


    //If DEM Error is desired - ensure that bperp and inc datasets work
    if (opts->estimateDEMError)
    {
        incDataset = reinterpret_cast<GDALDataset*>(GDALOpen( opts->incAngleFile.c_str(), GA_ReadOnly));
        if ((incDataset == NULL))
        {
            std::cout << "Cannot open file { " << opts->incAngleFile << "} for reading. \n";
            for (int jj=0; jj < opts->nPairs; jj++)
            {
                GDALClose( unwDataset[jj]);
            }
            GDALDestroyDriverManager();

            throw std::invalid_argument("Invalid incidence angle encountered");
        }

        if ((incDataset->GetRasterXSize() != inCols) || 
            (incDataset->GetRasterYSize() != inRows))
        {
            for(int jj=0; jj<opts->nPairs; jj++)
            {
                GDALClose(unwDataset[jj]);
            }
            GDALClose(incDataset);
            GDALDestroyDriverManager();
            throw std::invalid_argument("Inc file size mismatch");
        }


        for(int ii=0; ii<opts->nSAR; ii++)
        {
            bperpDataset[ii] = reinterpret_cast<GDALDataset*>(GDALOpen( scenes[ii].bperpName.c_str(), GA_ReadOnly));
            if (bperpDataset[ii] == NULL)
            {
                std::cout << "Cannot open file { " << scenes[ii].bperpName << "} for reading \n";
                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj<ii; jj++)
                {
                    GDALClose(bperpDataset[jj]);
                }
                GDALClose(incDataset);
                GDALDestroyDriverManager();

                throw std::invalid_argument("Invalid bperp file encountered");
            }

            if ((bperpDataset[ii]->GetRasterXSize() != inCols) || 
                (bperpDataset[ii]->GetRasterYSize() != inRows))
            {
                for(int jj=0; jj< opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }
                GDALClose(incDataset);

                for(int jj=0; jj<=ii; jj++)
                {
                    GDALClose(bperpDataset[jj]);
                }
                GDALDestroyDriverManager();
                throw std::invalid_argument("Bperp file size mismatch");
            }
        }
    } 

    //Compute output dimensions
    int nRows = opts->bbox[3];
    int nCols = opts->bbox[2];

    std::cout << "Output raster size: " << nCols << " P x " << nRows << " L\n";

    //Start setting up output directories and files
    std::vector< GDALDataset*> outDataset(opts->nSAR);
    GDALDataset* tcorrDataset;
    GDALDataset* velDataset;
    GDALDataset* zerrDataset;

    {
        //Output driver properties
        GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
        char **mOptions = NULL;
        mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
        mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

        {
            for(int ii=0; ii<opts->nSAR; ii++)
            {
                const char* fname = CPLStrdup(CPLFormFilename(opts->outDir.c_str(), opts->dates[ii].c_str(), ".bin"));
    
                outDataset[ii] = reinterpret_cast<GDALDataset*> (poDriver->Create(fname, nCols, nRows, 1, GDT_Float32, mOptions));

                if (outDataset[ii] == NULL)
                {
                    std::cout << "Could not create output file: " << fname << "\n";
               
                    for(int jj=0; jj< opts->nPairs; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }

                    for(int jj=0; jj<ii; jj++)
                    {
                        GDALClose(outDataset[jj]);
                    }

                    if (opts->estimateDEMError)
                    {
                        GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }

                    GDALDestroyDriverManager();

                    throw std::invalid_argument("Could not create output time series file");
                }      
            }

            const char* fname = CPLStrdup(CPLFormFilename(opts->outDir.c_str(), "tcorr",".bin"));
            tcorrDataset = reinterpret_cast<GDALDataset*> (poDriver->Create(fname, nCols, nRows, 1, GDT_Float32, mOptions));

            if (tcorrDataset == NULL)
            {
                std::cout << "Could not create output file: " << fname << "\n";

                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj< opts->nSAR; jj++)
                {
                    GDALClose(outDataset[jj]);
                }

                if (opts->estimateDEMError)
                {
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();

                throw std::invalid_argument("Could not create output tcorr file");
            }
        }


        {
            const char* fname = CPLStrdup(CPLFormFilename(opts->outDir.c_str(), "velocity",".bin"));
            velDataset = reinterpret_cast<GDALDataset*> (poDriver->Create(fname, nCols, nRows, 1, GDT_Float32, mOptions));

            if (velDataset == NULL)
            {
                std::cout << "Could not create output file: " << fname << "\n";

                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj< opts->nSAR; jj++)
                {
                    if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                }

                GDALClose(tcorrDataset); 
                
                if (opts->estimateDEMError)
                {
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();

                throw std::invalid_argument("Could not create output tcorr file");
            }
        }


        if (opts->estimateDEMError)
        {
            const char* fname = CPLStrdup(CPLFormFilename(opts->outDir.c_str(), "zerr",".bin"));
            zerrDataset = reinterpret_cast<GDALDataset*> (poDriver->Create(fname, nCols, nRows, 1, GDT_Float32, mOptions));

            if (zerrDataset == NULL)
            {
                std::cout << "Could not create output file: " << fname << "\n";

                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj< opts->nSAR; jj++)
                {
                    if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                }

                GDALClose(tcorrDataset);
                GDALClose(velDataset);

                if (opts->estimateDEMError)
                {
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();

                throw std::invalid_argument("Could not create output zerr file");
            }
        }
    }


    //Actual processing block 
    int blockysize = 128;
    if (blockysize > nRows)
        blockysize = nRows;

    std::cout << "Block size = " << blockysize << " lines \n";

    //Variables for handling blocks
    int blockcount = 0;
    int nBlocks = ((nRows-1)/blockysize) + 1;
   
    //Data variables
    arma::fmat obs(blockysize*nCols, opts->nPairs);
    arma::fmat rawts(blockysize*nCols, opts->nSAR-1);

    arma::fvec unwdata(blockysize*nCols);
    arma::fvec veldata(blockysize*nCols);
    arma::fvec tcorrdata(blockysize*nCols);


    arma::fmat bperp;
    arma::fvec bperpref;
    arma::fvec zerrdata;
    arma::fvec incdata;
    if (opts->estimateDEMError)
    {
        bperp.resize(blockysize*nCols, opts->nSAR-1);
        zerrdata.resize(blockysize*nCols);
        incdata.resize(blockysize*nCols);
        bperpref.resize(blockysize*nCols);
    }


    for(int yoff=0; yoff < nRows; yoff += blockysize)
    {
        //Increment block counter
        blockcount++;

        std::cout << "Processing block " << blockcount << " / " << nBlocks << "\n";

        double t_start = getWallTime();

        //Determine number of lines to read
        int inysize = blockysize;
        if ((yoff+inysize) > nRows)
            inysize = nRows - yoff;

        veldata.zeros();
        rawts.zeros();
        tcorrdata.zeros();
        obs.zeros();

        std::cout << "Reading unwrapped data \n"; 
        for(int ii=0; ii<opts->nPairs; ii++)
        {

            int status = unwDataset[ii]->GetRasterBand(1)->RasterIO(GF_Read,
                        opts->bbox[0], opts->bbox[1]+yoff,
                        nCols, inysize,
                        unwdata.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

            if (status != 0)
            {
                std::cout << "Error reading block at " << yoff 
                          << " for pair " << opts->pairs[ii].ifgName 
                          << "\n";

                for(int jj=0; jj< opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj< opts->nSAR; jj++)
                {
                    if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                }

                GDALClose(tcorrDataset);
                GDALClose(velDataset);

                if (opts->estimateDEMError)
                {
                    GDALClose(zerrDataset);
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();
                throw std::invalid_argument("Error reading unwrapped data");
            }

            obs.col(ii) = unwdata - opts->pairs[ii].referenceOffset;
            GDALTermProgress(ii / (1.0 * opts->nPairs), NULL, NULL);
        }

        double t_end = getWallTime();
        std::cout << "\n Time to read unwrapped data block: " << (t_end-t_start)/60.0 << "\n";

        //Reading incidence angle and bperp values
        if (opts->estimateDEMError)
        {

            t_start = getWallTime();
            zerrdata.zeros();
            incdata.zeros();
            bperpref.zeros();
            bperp.zeros();

            {
                //Read in the incidence angle
                int status = incDataset->GetRasterBand(1)->RasterIO(GF_Read,
                        opts->bbox[0], opts->bbox[1]+yoff,
                        nCols, inysize,
                        incdata.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

                if (status != 0)
                {
                    std::cout << "Error reading inc block at " << yoff 
                          << "\n";

                    for(int jj=0; jj< opts->nPairs; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }

                    for(int jj=0; jj< opts->nSAR; jj++)
                    {
                        if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                    }

                    GDALClose(tcorrDataset);
                    GDALClose(velDataset);

                    if (opts->estimateDEMError)
                    {
                        GDALClose(zerrDataset);
                        GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }

                    GDALDestroyDriverManager();
                    throw std::invalid_argument("Error reading inc angle data");
                }

                std::cout << "inc block: " << arma::min(incdata) << "  " << arma::max(incdata) << "\n";
                incdata = arma::sin((M_PI/180.0) * incdata);
            }
            
            {
                //First read in refIndex
                int status = bperpDataset[sbasInv.refIndex]->GetRasterBand(1)->RasterIO(GF_Read,
                        opts->bbox[0], opts->bbox[1]+yoff,
                        nCols, inysize,
                        bperpref.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

                if (status != 0)
                {
                    std::cout << "Error reading bperp block at " << yoff 
                          << " for date " << scenes[sbasInv.refIndex].date
                          << "\n";

                    for(int jj=0; jj< opts->nPairs; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }

                    for(int jj=0; jj< opts->nSAR; jj++)
                    {
                        if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                    }

                    GDALClose(tcorrDataset);
                    GDALClose(velDataset);

                    if (opts->estimateDEMError)
                    {
                        GDALClose(zerrDataset);
                        GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }

                    GDALDestroyDriverManager();
                    throw std::invalid_argument("Error reading bperp data");
                }
            }

            for(int ii=0; ii< opts->nSAR-1; ii++)
            {
                int ind = ii + (ii >= sbasInv.refIndex);

                int status = bperpDataset[ind]->GetRasterBand(1)->RasterIO(GF_Read,
                        opts->bbox[0], opts->bbox[1]+yoff,
                        nCols, inysize,
                        bperp.colptr(ii),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

                if (status != 0)
                {
                    std::cout << "Error reading bperp block at " << yoff 
                          << " for date " << scenes[ii].date
                          << "\n";

                    for(int jj=0; jj< opts->nPairs; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }

                    for(int jj=0; jj< opts->nSAR; jj++)
                    {
                        if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                    }

                    GDALClose(tcorrDataset);
                    GDALClose(velDataset);

                    if (opts->estimateDEMError)
                    {
                        GDALClose(zerrDataset);
                        GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }

                    GDALDestroyDriverManager();
                    throw std::invalid_argument("Error reading bperp data");
                }


                bperp.col(ii) = bperp.col(ii) - bperpref;

                std::cout << "Bperp block: " << ii << "  " << arma::min(bperp.col(ii)) << "  " << arma::max(bperp.col(ii)) << "\n";

                double dr = opts->rangeSpacing;
                double r0 = opts->startingRange + dr * opts->bbox[0];
                double wvl = opts->wavelength;

            #pragma omp parallel for\
                default(shared)
                for(int jj=0; jj < (inysize*nCols); jj++)
                {
                    double rng = r0 + (jj%nCols) * dr;
                    bperp.at(jj,ii) *= (4.0 * M_PI)/(wvl * rng * incdata[jj]);
                }


                GDALTermProgress((ii+2)/(opts->nSAR+1.0), NULL, NULL);
            }

            t_end = getWallTime();
            std::cout <<"\nTime to read in Bperp data: " << (t_end - t_start)/60.0 << "\n";
        }


        //Block for executing SBAS
        {
            t_start = getWallTime();

            Solution workers[nThreads];
            for(int jj=0; jj<nThreads; jj++)
            {
                workers[jj].tsest.resize(opts->nSAR-1);
            }

            for(int jj=0; jj<inysize; jj++)
            {

                bool getDEMError = opts->estimateDEMError;

        #pragma omp parallel for\
                default(shared)
                for(int kk=0; kk<nCols;kk++)
                {
                    int threadnum = omp_get_thread_num();

                    int ind = jj*nCols + kk;

                    //Extract row of data
                    arma::frowvec dph = obs.row(ind);
                    
                    //Make it look like a column vector
                    arma::fvec dphcol(dph.memptr(), dph.n_elem, false);

                    //Always Solve SBAS equation
                    sbasInv.solve(dphcol, &(workers[threadnum])); 

                   
                    if(getDEMError)
                    {
                        //Estimate velocity in combination with dem error
                        arma::frowvec bperprow = bperp.row(ind);
                        arma::fvec bperpcol(bperprow.memptr(), bperprow.n_elem, false);

                        sbasInv.estimateDEMError(bperpcol, &(workers[threadnum]));
                        zerrdata[ind] = workers[threadnum].zerr;
                    }
                    else
                    {
                        //Else just fit a model
                        sbasInv.estimateVelocity(&(workers[threadnum]));
                    }

                    veldata[ind] = workers[threadnum].velocity; 
                    tcorrdata[ind] = workers[threadnum].temporalCorrelation;

                    //Copy testrow back as a row
                    rawts.row(ind) =  arma::frowvec(workers[threadnum].tsest.memptr(), opts->nSAR-1, false); 

                }

                GDALTermProgress(jj / (1.0*inysize), NULL, NULL);
            }

            t_end = getWallTime();
            std::cout << "\n Time to process block: " << (t_end-t_start)/60.0 << "\n";
        }

        //Block for writing outputs
        if (opts->estimateDEMError)
        {
            int status = zerrDataset->GetRasterBand(1)->RasterIO(GF_Write,
                        0, yoff,
                        nCols, inysize,
                        zerrdata.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

            if (status != 0)
            {
                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj< opts->nSAR; jj++)
                {
                    if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                }

                GDALClose(tcorrDataset);
                GDALClose(velDataset);
                if (opts->estimateDEMError)
                {
                    GDALClose(zerrDataset);
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();
                throw std::invalid_argument("Error writing data to dem error file");
            }
        }

        {

            int status = velDataset->GetRasterBand(1)->RasterIO(GF_Write,
                        0, yoff,
                        nCols, inysize,
                        veldata.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

            if (status != 0)
            {
                for(int jj=0; jj<opts->nPairs; jj++)
                {
                    GDALClose(unwDataset[jj]);
                }

                for(int jj=0; jj < opts->nSAR; jj++)
                {
                    if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                }

                GDALClose(tcorrDataset);
                GDALClose(velDataset);

                if (opts->estimateDEMError)
                {
                    GDALClose(zerrDataset);
                    GDALClose(incDataset);
                    for(int jj=0; jj<opts->nSAR; jj++)
                    {
                        GDALClose(bperpDataset[jj]);
                    }
                }

                GDALDestroyDriverManager();
                throw std::invalid_argument("Error writing data to count file");
            }

            {
                for(int jj=0; jj<opts->nSAR-1;jj++)
                {
                    int ind = jj + (jj>=sbasInv.refIndex);
                    
                    int status = outDataset[ind]->GetRasterBand(1)->RasterIO(GF_Write,
                            0, yoff,
                            nCols, inysize,
                            rawts.colptr(jj),
                            nCols, inysize, GDT_Float32,
                            sizeof(float), sizeof(float)*nCols, NULL);

                    if (status != 0)
                    {
                        for(int kk=0; kk<opts->nPairs; kk++)
                        {
                            GDALClose(unwDataset[kk]);
                        }

                        for(int kk=0; kk < opts->nSAR; kk++)
                        {
                            if (outDataset[kk] != NULL) GDALClose(outDataset[kk]);
                        }

                        GDALClose(tcorrDataset);
                        GDALClose(velDataset);
                        GDALClose(zerrDataset);
                        if (opts->estimateDEMError)
                        {
                            GDALClose(zerrDataset);
                            GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }


                        GDALDestroyDriverManager();
                        throw std::invalid_argument("Error writing ts estimates");
                    }
                }

                int status = tcorrDataset->GetRasterBand(1)->RasterIO(GF_Write,
                        0, yoff,
                        nCols, inysize,
                        tcorrdata.memptr(),
                        nCols, inysize, GDT_Float32,
                        sizeof(float), sizeof(float)*nCols, NULL);

                if (status != 0)
                {
                    for(int jj=0; jj<opts->nPairs; jj++)
                    {
                        GDALClose(unwDataset[jj]);
                    }

                    for(int jj=0; jj < opts->nSAR; jj++)
                    {
                        if (outDataset[jj] != NULL) GDALClose(outDataset[jj]);
                    }

                    GDALClose(tcorrDataset);
                    GDALClose(velDataset);
                    if (opts->estimateDEMError)
                    {
                        GDALClose(zerrDataset);
                        GDALClose(incDataset);
                        for(int jj=0; jj<opts->nSAR; jj++)
                        {
                            GDALClose(bperpDataset[jj]);
                        }
                    }

                    GDALDestroyDriverManager();
                    throw std::invalid_argument("Error writing data to count file");
                }

            }
        }
    }


    //Close all input datasets
    for(int ii=0; ii<opts->nPairs; ii++)
    {
        GDALClose(unwDataset[ii]);
    }

    if (opts->estimateDEMError)
    {
        GDALClose(zerrDataset);
        GDALClose(incDataset);
        for(int ii=0; ii<opts->nSAR; ii++)
        {
            GDALClose(bperpDataset[ii]);
        }
    }

    for(int ii=0; ii<opts->nSAR; ii++)
    {
        GDALClose(outDataset[ii]);
    }

    GDALClose(tcorrDataset);
    GDALClose(velDataset);

    GDALDestroyDriverManager();
}

