/**\file ampdispersion.c
 * \author Piyush Agram.
 *  */

#include "evd.hpp"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "cpl_port.h"
#include "fringe/fringe_common.hpp"
#include "fringe/EigenLapack.hpp"

#include <vector>
#include <armadillo>

int evd_process(evdOptions *opts)
{

    //Print user options to screen
    opts->print();


    int cols, rows, nbands;
    int blockysize, boxesperblock;
    int miniStackCount = opts->miniStackCount;
    
    //First thing to do is to read dimenions.
    GDALDataset* inDataset = NULL;          //Input stack
    GDALDataset* wtsDataset = NULL;         //Input weights
    GDALDataset* corrDataset = NULL;        //Output correlation dataset
    GDALDataset* compDataset = NULL;        //Compressed SLC dataset
    std::vector<GDALDataset*> outDataset;   //Vector of output SLCs
    
    //Vector of dates for internal bookkeeping
    std::vector<std::string> dates;

    //Clock variables
    double t_start, t_end;

    //Make sure C++11 is being used
    //unsigned int - 32 bytes
    checkLongSetting();

    //Register GDAL drivers
    GDALAllRegister();

    //Determine the number of uint32 bytes needed to store weights
    int numAround = (2*opts->Ny + 1) * (2*opts->Nx + 1); 
    int nulong = ceil( (numAround*1.0)/ 32.0);
    std::cout << "Number of uint32 bytes for mask: " << nulong << "\n";

    //Open the stack dataset
    inDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->inputDS.c_str(), GA_ReadOnly));
    if (inDataset == NULL)
    {
        std::cout << "Cannot open stack file { " <<  opts->inputDS << " } with GDAL for reading. \n";

        std::cout << "GDALOpen failed - " << opts->inputDS << "\n";
        std::cout << "Exiting with error code .... (102) \n";
        GDALDestroyDriverManager();
        return 102;
    }

    cols = inDataset->GetRasterXSize();
    rows = inDataset->GetRasterYSize();
    nbands = inDataset->GetRasterCount();

    std::cout << "Number of rows  = " << rows << "\n";
    std::cout << "Number of cols  = " << cols << "\n";
    std::cout << "Number of bands = " << nbands << "\n";


    //Resolve method and bandwidth
    if (opts->method.compare("STBAS") == 0)
    {
        //If bandwidth is not provided
        if (opts->bandWidth <= 0)
        {
            std::cout << "Requested STBAS but no bandwidth provided \n";
            GDALDestroyDriverManager();
            return 101;
        }
        else if (opts->bandWidth >= (nbands-1))
        {
            std::cout <<"Requested STBAS bandwidth " 
                      << opts->bandWidth << " is larger than full bandwidth " 
                      << (nbands-1) << "\n";
            GDALDestroyDriverManager();
            return 101;
        }
    }

    //Open weights file
    wtsDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->wtsDS.c_str(), GA_ReadOnly));
    if (wtsDataset == NULL)
    {
        std::cout << "Could not open weights dataset {" << opts->wtsDS << "}\n";
        std::cout << "Exiting with non-zero error code ... 105 \n";
        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 105;
    }

    //Check for consistency between input and weight file
    {
        int code = 0;
        int cc = wtsDataset->GetRasterXSize();
        if (cc != cols)
        {
            std::cout << "Width mismatch between input dataset and weight dataset\n";
            code = 106;
        }

        int rr = wtsDataset->GetRasterYSize();
        if (rr != rows)
        {
            std::cout << "Length mismatch between input dataset and weight dataset \n";
            code = 107;
        }

        int bb = wtsDataset->GetRasterCount();
        if (bb != nulong)
        {
            std::cout << "Number of bands mismatch for weights and window size \n";
            code = 108;
        }

        {
            int inNx = ::atoi(wtsDataset->GetMetadataItem("HALFWINDOWX", "ENVI"));
            int inNy = ::atoi(wtsDataset->GetMetadataItem("HALFWINDOWY", "ENVI"));

            if (inNx != opts->Nx) 
            {
                std::cout << "Half window size x of wts is different from input. \n";
                code = 109;
            }

            if (inNx == 0)
            {
                std::cout << "No non-zero metadata item called HALFWINDOWX \n";
                code = 109;
            }

            if (inNy != opts->Ny)
            {
                std::cout << "Half window size y of wts is different from input. \n";
                code = 110;
            }

            if (inNy == 0)
            {
                std::cout << "No non-zero metadata item called HALFWINDOWY \n";
                code = 110;
            }
        }

        if (code != 0)
        {
            std::cout <<"Exiting with error code ....("<<code<<")\n";
            GDALClose(inDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();

            return code;
        }
    }


    //Determine blocksizes
    boxesperblock = int((opts->memsize * 1.0e6)/cols) / (opts->blocksize * (nbands*20+4+nulong)) ;
    blockysize = boxesperblock * opts->blocksize;
    if (blockysize < opts->blocksize)
    {
        blockysize = opts->blocksize;
        boxesperblock = 1;
    }

    if (blockysize > rows)
    {
        blockysize = rows;
        boxesperblock = 1;
    }

    std::cout << "Block size = " << blockysize << " lines \n";

    int totalblocks = ceil( rows / (1.0 * blockysize));
    std::cout << "Total number of blocks to process: " << totalblocks << "\n";

    //Start the clock
    t_start = getWallTime();

    //Temporary array for reading in complex data
    arma::cx_fmat cpxdata(cols*blockysize, nbands); //To store input data
    arma::cx_fmat evddata(cols*blockysize, nbands);  //To store output data
    arma::Mat<unsigned int> wts(nulong, cols*blockysize); //To store weights
    arma::fvec corr(cols*blockysize);    //To store correlation
    arma::cx_fvec comp(cols*blockysize); //To store compressed SLC

    //Start block-by-block processing
    int blockcount = 0;
    int nthreads = 0;

    //Resize vectors based on number of bands
    outDataset.resize(nbands);
    dates.resize(nbands);
  
    //Creation of output datasets
    {
        //Make sure that all the input stack layers have dates
        for(int b=1; b<=nbands; b++)
        {
            const char* date = inDataset->GetRasterBand(b)->GetMetadataItem("Date", "slc");
            if (strnlen(date, 9) != 8)
            {
                std::cout << "Band " << b << " does not appear to have Date information in slc metadata domain \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALDestroyDriverManager();

                return 112;
            }
            dates[b-1] =  std::string(date);
        }

        //Check output folder
        VSIStatBufL sStat;

        //If folder already exists
        if ( VSIStatExL( opts->outputFolder.c_str(), &sStat, VSI_STAT_EXISTS_FLAG | VSI_STAT_NATURE_FLAG) == 0)
        {
            if ( VSI_ISDIR(sStat.st_mode) )
            {
                std::cout << "Output folder : " << opts->outputFolder << " already exists \n";
                std::cout << "Returning without processing. Clean up output folder and rerun.. \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALDestroyDriverManager();

                return 113;
            }
            else
            {
                std::cout << opts->outputFolder << " already exists and appears to be a file on disk. \n";
                std::cout << "Returning without processing. Clean up output folder and rerun.. \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALDestroyDriverManager();

                return 114;
            }
        }

        //Create output folder
        if (VSIMkdir( opts->outputFolder.c_str(), 0777) != 0)
        {
            std::cout << "Could not create output folder: " << opts->outputFolder << "\n";

            GDALClose(inDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();
            return 115;
        }
        std::cout << "Created output folder: " << opts->outputFolder << "\n";

        //Creation of output dataset 
        for(int b=1; b<=nbands; b++)
        {
            GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI"); 
            char **mOptions = NULL;
            mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
            mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

            std::string fname = CPLFormFilename(opts->outputFolder.c_str(), dates[b-1].c_str(), ".slc");
            std::cout << "Input: " << fname <<  " " << fname.size() <<  "\n";

            GDALDataset *temp = (GDALDataset*) poDriver->Create(fname.c_str(), cols, rows, 1, GDT_CFloat32, mOptions);

            if (temp == NULL)
            {
                std::cout << "Could not create output SLC: " << fname << "\n";
                std::cout << "Exiting with non-zero error code ... 116 \n";

                GDALClose(inDataset);
                GDALClose(wtsDataset);

                for (int j=0; j < b-1; j++)
                {
                    GDALClose(outDataset[j]);
                }
                GDALDestroyDriverManager();

                return 116;
            }

            outDataset[b-1] = temp;

            std::cout << "Created : " << fname  <<  "  " << temp << "\n";
            
            CSLDestroy(mOptions);
        }


        //Temporal coherence
        {
            GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI"); 
            char **mOptions = NULL;
            mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
            mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

            std::string fname = CPLFormFilename(opts->outputFolder.c_str(), "tcorr.bin", NULL);

            std::cout << "Corr: " << fname <<  " " << fname.size() << "\n";
            corrDataset = (GDALDataset*) poDriver->Create(fname.c_str(), cols, rows, 1, GDT_Float32, mOptions);

            if (corrDataset == NULL)
            {
                std::cout << "Could not create temporal correlation file: " << fname << "\n";
                std::cout << "Exiting with non-zero error code ... 117 \n";

                GDALClose(inDataset);
                GDALClose(wtsDataset);

                for (int j=0; j<nbands; j++)
                {
                    GDALClose(outDataset[j]);
                }
                GDALDestroyDriverManager();
                return 117;
            }

            std::cout << "Created : " << fname << " " << corrDataset << "\n";

            CSLDestroy(mOptions);
        }

        //Compressed SLC
        {
            GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
            char **mOptions = NULL;
            mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
            mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

            std::string fname = CPLFormFilename(opts->outputCompressedSlcFolder.c_str(), opts->compSlc.c_str(), NULL);

            std::cout << "Compressed SLC: " << fname <<  " " << fname.size() << "\n";
            compDataset = (GDALDataset*) poDriver->Create(fname.c_str(), cols, rows, 1, GDT_CFloat32, mOptions);

            if (compDataset == NULL)
            {
                std::cout << "Could not create compressed SLC file: " << fname << "\n";
                std::cout << "Exiting with non-zero error code ... 117 \n";

                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(corrDataset);
                for (int j=0; j<nbands; j++)
                {
                    GDALClose(outDataset[j]);
                }
                GDALDestroyDriverManager();
                return 117;
            }

            std::cout << "Created : " << fname << " " << compDataset << "\n";
            CSLDestroy(mOptions);
        }



    }

    nthreads = numberOfThreads();

    //Useful variables
    int Ny = opts->Ny;
    int Nx = opts->Nx;
    int Wy = 2*Ny+1;
    int Wx = 2*Nx+1;

    const Ulongmask bitmask = {Ny,Nx};
    
    //Setup the Eigen Value workers
    EVWorker *workers = new EVWorker[nthreads];
    for(int ii=0; ii < nthreads; ii++)
    {
        workers[ii].prepare(nbands);
    }

    //Setup arrays for computing coherence matrix
    arma::cx_cube Numer(nbands, nbands, nthreads);
    arma::cx_cube Amp1(nbands, nbands, nthreads);
    arma::cx_cube Amp2(nbands,nbands, nthreads);
    arma::cx_cube Covar(nbands, nbands, nthreads);


    //Block-by-block processing
    int yoff = 0;

    while( yoff < rows)
    {
        //Increment block counter
        blockcount++;


        int inysize = blockysize;

        //Number of lines to read for the block
        if ((yoff+inysize) > rows)
            inysize = rows - yoff;

//        std::cout << "Block " << blockcount << ": " << inysize << " lines \n";

        //Block logic
        int firstlinetowrite;
        int linestowrite;
        int rollback;

        //For first block include the top half window.
        if (blockcount == 1)
        {
            firstlinetowrite = 0;
            linestowrite = inysize - Ny;
            rollback = Ny;
        
            if( (yoff+blockysize) >= rows) //There is only one block, write the whole thing
            {
                linestowrite = inysize;
                rollback = 0;
            }

        }
        else if( (yoff+blockysize) >= rows) //Last block. Write bottom half window.
        {
            firstlinetowrite = Ny;
            linestowrite = inysize-Ny;
            rollback = 0;
        }
        else
        {
            firstlinetowrite = Ny;
            linestowrite = inysize-2*Ny;
            rollback = 0;
        }


        //Init to zeros for each block
        cpxdata.zeros();
        evddata.zeros();
        wts.zeros();
        Numer.zeros();
        Amp1.zeros();
        Amp2.zeros();
        corr.zeros();
        comp.zeros();
        
        //Read in data from all SLCs
        for(int bb=0; bb < nbands; bb++)
        {
            int status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read,
                        0, yoff,   
                        cols, inysize,
                        (void*) (cpxdata.colptr(bb)), 
                        cols, inysize, GDT_CFloat32,
                        sizeof(std::complex<float>),
                        sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error reading data from band " << bb+1 << " at line " 
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (118) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                for(int jj=0; jj<nbands; jj++) GDALClose(outDataset[jj]);
                GDALClose(corrDataset);
                GDALClose(compDataset);
                GDALDestroyDriverManager();
                return 118;
            }
        }

        {
            //Read in the weights dataset
            int status =  wtsDataset->RasterIO(GF_Read,
                    0, yoff,
                    cols, inysize,
                    (void*) (wts.memptr()),
                    cols, inysize, GDT_UInt32,
                    nulong, NULL,
                    4*nulong, 4*nulong*cols,
                    4, NULL);
            if (status != 0)
            {
                std::cout << "Error reading weights at line " << yoff << "\n";
                std::cout << "Exiting with error code .... (119) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(corrDataset);
                GDALClose(compDataset);
                for (int jj=0; jj<nbands; jj++) GDALClose(outDataset[jj]);
                GDALDestroyDriverManager();
                return 119;
            }
        }

        //Copy method to constant string
        const std::string method = opts->method;
        const int BW = opts->bandWidth;
        const bool isstbas = (method.compare("STBAS") == 0);
        const bool ismle = (method.compare("MLE") == 0);

    #pragma omp parallel for\
        default(shared)
        for(int pp=(firstlinetowrite*cols); pp < ((firstlinetowrite+linestowrite)*cols); pp++)
        {
            unsigned int *cenptr = wts.colptr(pp);
            if (bitmask.getbit(cenptr,0,0) == 0) continue;

            //Get thread number
            int threadnum = omp_get_thread_num();

            Covar.slice(threadnum).zeros();
            Numer.slice(threadnum).zeros();
            Amp1.slice(threadnum).zeros();
            Amp2.slice(threadnum).zeros();

            int cenii = pp/cols;
            int cenjj = pp%cols;

            int xmin = std::max(cenjj-Nx, 0);
            int xmax = std::min(cols-1, cenjj+Nx);

            int ymin = std::max(cenii-Ny, 0);
            int ymax = std::min(inysize-1, cenii+Ny);


            //Gather data
            int npix = 0;
            for (int ii=ymin; ii <=ymax; ii++)
            {
                for (int jj=xmin; jj <= xmax; jj++)
                {
                    /*std::cout << "ii = " << ii - cenii
                              << ", jj = " << jj - cenjj
                              << ", flag = " << bitmask.getbit(cenptr, ii-cenii, jj-cenjj) << "\n";*/

                    if (bitmask.getbit(cenptr, ii-cenii, jj-cenjj) != 0)
                    {
                        npix++;
                        arma::cx_frowvec ptr = cpxdata.row(ii*cols + jj);

                        //Compute correlation for top half of the matrix only
                        for(int ti=0; ti < nbands; ti++)
                        {
                            for(int tj=ti+1; tj<nbands; tj++)
                            {
                                Numer.at(ti,tj,threadnum) += ptr[ti] * std::conj(ptr[tj]);
                                Amp1.at(ti,tj,threadnum)  += std::complex<double>(std::pow(std::abs(ptr[ti]),2),0.);
                                Amp2.at(ti,tj,threadnum)  += std::complex<double>(std::pow(std::abs(ptr[tj]),2),0.);
                            }
                        }
                    }
                }
            }

            if (npix < 2) continue;     //Minimum number of neighbors

            //Create coherence matrix - common to all methods
            for (int ti=0; ti <nbands; ti++)
            {
                for(int tj=ti+1; tj< nbands; tj++)
                {
                    std::complex<double> val1 = Numer.at(ti,tj,threadnum);
                    double a1 = Amp1.at(ti,tj,threadnum).real();
                    double a2 = Amp2.at(ti,tj,threadnum).real();

                    std::complex<double> res = val1 / std::sqrt(a1*a2);
                    Covar.at(ti,tj,threadnum) = res;
                    Covar.at(tj,ti,threadnum) = std::conj(res);
                }
                Covar.at(ti,ti,threadnum) = std::complex<double>(1.0, 0.0);
            }

            if (ismle)
            {
                //Store absolute value in Amp1
                for (int ti=0; ti<nbands; ti++)
                {
                    for(int tj=ti+1; tj<nbands; tj++)
                    {
                        double coh = std::abs( Covar.at(ti,tj,threadnum));
                        Amp1.at(ti,tj,threadnum) = std::complex<double>(coh, 0.0);
                        Amp1.at(tj,ti,threadnum) = std::complex<double>(coh, 0.0);
                    }
                    Amp1.at(ti,ti,threadnum) = std::complex<double>(1.0, 0.0);
                }

                bool valid = true;
                //This is check to see if coherence is numerically well behaved
                {
                    //Copy covar to amp2
                    Amp2.slice(threadnum) = Covar.slice(threadnum);

                    //Check if covariance positive semi definite
                    int status = workers[threadnum].smallestEigen(Amp2.slice_memptr(threadnum), false);
                    if (status != 0)
                    {
                        corr[pp] = -1;
                        valid = false;
                    }
                    else
                    {
                        if (workers[threadnum].eigval[0] < 1.0e-6)
                        {
                            corr[pp] = -2;
                            valid = false;
                        }
                    }
                }
                if (!valid) continue;   //Skip rest of computation

                //This is check to see if coherence is numerically well behaved
                {
                    //Copy coherence to Amp2
                    Amp2.slice(threadnum) = Amp1.slice(threadnum);

                    //Check if coherence is positive semi-definite
                    int status = workers[threadnum].smallestEigen(Amp2.slice_memptr(threadnum), false);

                    if (status !=0) 
                    {
                        corr[pp] = -3;
                        valid = false;
                    }
                    else
                    {
                        if (workers[threadnum].eigval[0] < 1.0e-6)
                        {
                            corr[pp] = -4;
                            valid = false;
                        }
                    }
                }
                if (!valid) continue;   //Skip rest of computation

                //MLE stuff
                //We have access to eigenvector as well as Covar here
                //This is where the iterative MLE solver goes
                //MLE solver should put solution back into evddata

                //Compute inverse of the coherence matrix here
                {
                    //copy coherence to Amp2
                    Amp2.slice(threadnum) = Amp1.slice(threadnum);

                    //Perform inverse operation
                    int status = workers[threadnum].positiveDefiniteInverse(Amp2.slice_memptr(threadnum));
                    if (status != 0)
                    {
                        corr[pp] = -5;
                        valid = false;
                    }
                }
                if (!valid) continue;   //Skip rest of computation

                //Perform Hadamard product
                Amp2.slice(threadnum) %= Covar.slice(threadnum);

                //Obtain the smallest eigen vector
                {
                    int status = workers[threadnum].smallestEigen( Amp2.slice_memptr(threadnum), true);
                    if (status != 0)
                    {
                        corr[pp] = -6;
                        valid = false;
                    }
                    else
                    {
                        if (workers[threadnum].eigval[0] < 1.0e-6)
                        {
                            corr[pp] = -7;
                            valid = false;
                        }
                    }
                }
                if (!valid) continue;   //Skip updating evddata

            } //End of MLE
            else
            {
                //This is for STBAS / EVD - which both use the same concept
                //Largest eigen vector of the correlation matrix
               
                //If STBAS, apply bandwidth limits
                if (isstbas)
                {
                    //Blank out covariance entries outside the bandwidth
                    for (int ti=0; ti <nbands; ti++)
                    {
                        for(int tj=ti+BW+1; tj< nbands; tj++)
                        {
                            Covar.at(tj,ti,threadnum) = std::complex<double>(0., 0.);
                            Covar.at(ti,tj,threadnum) = std::complex<double>(0., 0.);
                        }
                    }
                }

                //Copy covar into amp1 for computation
                //Need to preserve covar for temporal coherence computation
                Amp1.slice(threadnum) = Covar.slice(threadnum);

                //Obtain the largest eigen vector
                bool valid = true;
                {
                    int status = workers[threadnum].largestEigen( Amp1.slice_memptr(threadnum), true);
                    if (status != 0)
                    {
                        corr[pp] = -6;
                        valid = false;
                    }
                    else
                    {
                        if (workers[threadnum].eigval[0] < 1.0e-6)
                        {
                            corr[pp] = -7;
                            valid = false;
                        }
                    }
                }
                if (!valid) continue;   //Skip updating evddata
           
            } //End of EVD or STBAS


            //If we have gotten this far, "eigvec" should contain the solution
            //For MLE, its from smallestEigen and its largetstEigen for others
            //Ensure first date of acquired SLCs is set to zero
            {
                std::complex<double> cJ(0.0, 1.0);
                double ph0 = std::arg(workers[threadnum].eigvec[miniStackCount-1]);
                for(int ii=0; ii < nbands; ii++)
                {
                    double res = std::arg(workers[threadnum].eigvec[ii]);
                    evddata.at(pp, ii) = std::exp(cJ * (res - ph0));
                }

                //Ensure that the reference index is exactly 1 - i.e, phase 0
                evddata.at(pp, miniStackCount-1) = std::complex<float>(1.0, 0.0);
            }


            //compress acquired SLCs of the stack.
            //This is common to all methods
            //For sequential approach previous compressed SLCs are excluded in the current compression
            {
                std::complex<double> sumcomp = 0.0;
                for(int ii=miniStackCount-1; ii < nbands; ii++)
                {
                    sumcomp += cpxdata.at(pp, ii) * std::conj( evddata.at(pp,ii));
                }
                comp.at(pp) = sumcomp / (1.0 * (nbands - miniStackCount + 1));
            }            
            
            //This part might need some updating to iterate starting with the eigen value solution
            //Numerical simulations dont seem to show any impact - but all papers suggest this


            //Temporal coherence estimation
            //This is common to all methods
            {
                std::complex<double> tempcorr(0.0,0.0);
                std::complex<double> cJ(0.0,1.0);
                int counter = 0;
                for(int ti=0; ti<nbands; ti++)
                {
                    int ulim = isstbas ? (ti+BW+1) : nbands;
                    for(int tj=ti+1; tj<ulim; tj++)
                    {
                        int ind = ti + nbands * tj;
                    
                        tempcorr += std::exp(cJ * (std::arg(Covar.at(ti,tj,threadnum)) - std::arg( evddata.at(pp,ti)) + std::arg( evddata.at(pp,tj))));
                        counter++;
                    }
                }
                corr.at(pp) = std::abs(tempcorr) /(counter * 1.0);
            }

        }

       
        //Write the data to output bands
        for(int ii=0; ii<nbands; ii++)
        {
            int status = outDataset[ii]->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff+firstlinetowrite,
                        cols, linestowrite,
                        (void*) (evddata.colptr(ii) + firstlinetowrite*cols),
                        cols, linestowrite, GDT_CFloat32,
                        sizeof(std::complex<float>), sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error writing EVD data at line " << yoff << " for band " << ii+1 << "\n";
                std::cout << "Exiting with error code .... (120) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(corrDataset);
                GDALClose(compDataset);
                for(int jj=0; jj < nbands; jj++) GDALClose(outDataset[jj]);
                GDALDestroyDriverManager();
                return 120;
            }
        }

        //Write the temporal correlation data
        {
            int status = corrDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff+firstlinetowrite,
                    cols, linestowrite,
                    (void*) (corr.memptr()+firstlinetowrite*cols),
                    cols, linestowrite, GDT_Float32,
                    sizeof(float), sizeof(float)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error writing temporal correlation data at line " << yoff << "\n";
                std::cout << "Exiting with error code .... (121) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(corrDataset);
                GDALClose(compDataset);
                for(int jj=0; jj<nbands; jj++) GDALClose(outDataset[jj]);
                GDALDestroyDriverManager();
                return 121;
            }
        }


        //Write the compressed SLC data
        {
            int status = compDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff+firstlinetowrite,
                    cols, linestowrite,
                    (void*) (comp.memptr()+firstlinetowrite*cols),
                    cols, linestowrite, GDT_CFloat32,
                    sizeof(std::complex<float>), sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error writing temporal correlation data at line " << yoff << "\n";
                std::cout << "Exiting with error code .... (121) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(corrDataset);
                GDALClose(compDataset);
                for(int jj=0; jj<nbands; jj++) GDALClose(outDataset[jj]);
                GDALDestroyDriverManager();
                return 121;
            }
        }            


        //Update the counter
        if ((yoff+blockysize) < rows)
        {
            yoff += linestowrite - rollback;
        }
        else
        {
            yoff=rows;
        }

        GDALTermProgress( yoff/(1.0*rows), NULL, NULL);

    }

    //Free the workers
    delete [] workers;

    //Close the datasets
    GDALClose(inDataset);
    GDALClose(wtsDataset);
    for (int b=0; b<nbands;b++)
    {
        GDALClose(outDataset[b]);
    }
    GDALClose(corrDataset);
    GDALClose(compDataset);
    return(0);
       
};



/* Main driver*/
int main(int  argc, const char *argv[] ) {

    //Options  
    evdOptions opts;
    int status;

    //Parse command line options
    status = opts.initFromCmdLine(argc, argv);

    if (status != 0)
    {
        std::cout << "Error parsing command line inputs \n";
        std::cout << "Exiting with non-zero return code .... \n";
        return (101);
    }


    //Execute 
    status = evd_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing calamp \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);
       
};
