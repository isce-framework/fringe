/**\file ampdispersion.c
 * \author Piyush Agram.
 *  */

#include "despeck.hpp"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "fringe/fringe_common.hpp"
#include <armadillo>

int despeck_process(despeckOptions *opts)
{

    //Print user options to screen
    opts->print();

    int cols, rows, nbands;
    int blockysize, boxesperblock;


    //First thing to do is to read dimenions.
    GDALDataset* inDataset = NULL;
    GDALDataset* wtsDataset = NULL;
    GDALDataset* outDataset = NULL;

    //Clock variables
    double t_start, t_end;

    //Register GDAL drivers
    GDALAllRegister();

    //Output data type
    GDALDataType outType;

    //Determine the number of uint32 bytes needed to store weights
    int nulong = ceil( ((2*opts->Ny+1)*(2*opts->Nx+1))/ 32.0);
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

    if ((opts->ibands[0] <= 0) || (opts->ibands[0] > nbands))
    {
        std::cout << "Master band " << opts->ibands[0] << " outside the range of permissible bands\n";
        std::cout << "Exiting with error code ... (102) \n";
        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 102;
    }

    if ((opts->ibands[1] > 0) && (opts->ibands[1] > nbands))
    {
        std::cout << "Slave band " << opts->ibands[1] << "outside the range of permissible bands\n";
        std::cout << "Exiting with error code ... (102) \n";
        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 102;
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
    boxesperblock = int((opts->memsize * 1.0e6)/cols) / (opts->blocksize * 4 * (6+nulong)) ;
    blockysize = boxesperblock * opts->blocksize;
    if (blockysize < opts->blocksize)
    {
        blockysize = opts->blocksize;
        boxesperblock = 1;
    }

    std::cout << "Block size = " << blockysize << " lines \n";

    int totalblocks = ceil( rows / (1.0 * blockysize));
    std::cout << "Total number of blocks to process: " << totalblocks << "\n";

    //Start the clock
    t_start = getWallTime();

    //Temporary array for reading in complex data
    arma::cx_fvec cpxdata1(cols*blockysize);
    arma::cx_fvec cpxdata2(cols*blockysize);
    arma::Mat<unsigned int> wts(nulong, cols*blockysize);
    arma::cx_fvec filt(cols*blockysize);

    //Start block-by-block processing
    int blockcount = 0;
    int status;
    int nthreads;


    //Creation of output dataset 
    {
        GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
        char **mOptions = NULL;
        mOptions = CSLSetNameValue( mOptions, "INTERLEAVE", "BIP");
        mOptions = CSLSetNameValue( mOptions, "SUFFIX", "ADD");
    
        outType = GDT_Float32;
        if (opts->ibands[1] > 0)
        {
            outType = GDT_CFloat32;
        }
        
        outDataset = (GDALDataset*) poDriver->Create(opts->outputDS.c_str(), cols, rows, 1, outType, mOptions); 
        
        if (outDataset == NULL)
        {
            std::cout << "Could not create despecked dataset {" << opts->outputDS << "} \n";
            std::cout << "Exiting with non-zero error code ... 104 \n";

            GDALClose(inDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();
            return 104;
        }
        CSLDestroy(mOptions);
    }
    

    nthreads = numberOfThreads();

    //Useful variables
    int Ny = opts->Ny;
    int Nx = opts->Nx;
    int Wy = 2*Ny+1;
    int Wx = 2*Nx+1;
    int band1 = opts->ibands[0];
    int band2 = opts->ibands[1];



    const Ulongmask bitmask = {Ny,Nx};

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

        //Init to zeros
        cpxdata1.zeros();
        cpxdata2.zeros();
        wts.zeros();
        
        
        //Read in the data
        status = inDataset->GetRasterBand(band1)->RasterIO( GF_Read,
                    0, yoff,   
                    cols, inysize,
                    (void*) (cpxdata1.memptr()), 
                    cols, inysize, GDT_CFloat32,
                    sizeof(std::complex<float>),
                    sizeof(std::complex<float>)*cols, NULL);

        if (status != 0)
        {
            std::cout << "Error reading data from band " << band1 << " at line " 
                    << yoff << "\n";
            std::cout << "Exiting with error code .... (108) \n";
            GDALClose(inDataset);
            GDALClose(wtsDataset);
            GDALClose(outDataset);
            GDALDestroyDriverManager();
            return 108;
        }


        if (band2 > 0)
        {
            status = inDataset->GetRasterBand(band2)->RasterIO( GF_Read,
                        0, yoff,
                        cols, inysize,
                        (void*) (cpxdata2.memptr()),
                        cols, inysize, GDT_CFloat32,
                        sizeof(std::complex<float>),
                        sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {   
                std::cout << "Error reading data from band " << band2 << " at line " 
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(inDataset);
                GDALClose(wtsDataset);
                GDALClose(outDataset);
                GDALDestroyDriverManager();
                return 108;
            }

            //If coherence is desired instead of amplitude
            if (opts->computeCoherence)
            {

    #pragma omp parallel for\
        default(shared)
                for(int jj=0; jj< (inysize*cols); jj++)
                {
                    float amp1 = std::abs(cpxdata1(jj));
                    float amp2 = std::abs(cpxdata2(jj));
                    cpxdata1(jj) *= std::conj( cpxdata2(jj));
                    cpxdata2(jj) = std::complex<float>(amp1*amp1, amp2*amp2);
                }
            }
            else
            {

    #pragma omp parallel for\
        default(shared)
                for(int jj=0; jj< (inysize*cols); jj++)
                {
                    cpxdata1(jj) *= std::conj(cpxdata2(jj));
                    cpxdata2(jj) = 1.0;
                }
            }
        }
        else
        {

        //If only a single band is provided, changing the data to abs values
    #pragma omp parallel for\
            default(shared)
            for(int jj=0; jj< (inysize*cols); jj++)
            {
                cpxdata1(jj) = std::abs(cpxdata1(jj));
            }
            cpxdata2.ones();
        }

        filt.zeros();

        //Read in the wts
        status = wtsDataset->RasterIO( GF_Read,
                    0, yoff,
                    cols, inysize,
                    (void*) (wts.memptr()),
                    cols, inysize, GDT_UInt32,
                    nulong, NULL,
                    4 * nulong, 4 * nulong * cols,
                    4, NULL);

        if (status != 0)
        {
            std::cout << "Error reading weights band at line "<< yoff << "\n";
            std::cout << "Exiting with error code .... (108) \n";
            GDALClose(inDataset);
            GDALClose(outDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();
            return 108;
        }


    #pragma omp parallel for\
        default(shared)
        for(int pp=(firstlinetowrite*cols); pp< ((firstlinetowrite+linestowrite)*cols); pp++)
        {
            unsigned int *cenptr = wts.colptr(pp);
            if (bitmask.getbit(cenptr, 0, 0) == 0) continue;

            int cenii = pp/cols;
            int cenjj = pp%cols;

            int xmin = std::max(cenjj-Nx, 0);
            int xmax = std::min(cols-1, cenjj+Nx);

            int ymin = std::max(cenii-Ny, 0);
            int ymax = std::min(inysize-1, cenii+Ny);

            std::complex<float> val = 0.0;
            std::complex<float> sumw = 0.0;

            for(int ii=ymin; ii<=ymax; ii++)
            {
                for(int jj=xmin; jj<=xmax; jj++)
                {
                    if (bitmask.getbit(cenptr, ii-cenii, jj-cenjj)!=0)
                    {
                        val += cpxdata1(ii*cols + jj);
                        sumw += cpxdata2(ii*cols + jj);
                    }
                }
            }

            //Summing and normalization
            if (sumw.real() > 0)
            {
                if (opts->computeCoherence)
                {
                    if (sumw.imag() > 0)
                    {
                        filt(pp) = val / (std::sqrt(sumw.real()) * std::sqrt(sumw.imag()));
                    }
                }
                else
                {
                    filt(pp) = val / sumw.real();
                }
            }
        }



        //Print debugs
/*        std::cout << "Block number: " << blockcount << "\n";
        std::cout << "yoff: " << yoff << "\n";
        std::cout << "firstlinetowrite: " << firstlinetowrite << "\n";
        std::cout << "linestowrite: " << linestowrite << "\n";
        std::cout << "rollback: " << rollback << "\n";

        std::cout << "Gdal block \n";
        std::cout << "First line : " << yoff + firstlinetowrite << "\n";
        std::cout << "Number lines: " << linestowrite << "\n";
        std::cout << "Next read start: " << yoff + linestowrite - rollback << "\n"; */

        //Write stuff to output files


        //Write to count file
        status = outDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff+firstlinetowrite,
                    cols, linestowrite,
                    (void*) (&filt[firstlinetowrite*cols]), 
                    cols, linestowrite, outType,
                    sizeof(std::complex<float>), sizeof(std::complex<float>)*cols, NULL);

        if (status != 0)
        {
            std::cout << "Error writing despeck data at line " << yoff << "\n";
            std::cout << "Exiting with error code .... (110) \n";
            GDALClose(inDataset);
            GDALClose(outDataset);
            GDALClose(wtsDataset);
            GDALDestroyDriverManager();
            return 110;
        }


//        std::cout << "Weights at (10,10): " << count(10,10) << " " 
//                  << wts(3,10*cols+10)<< "\n";
//        bitmask.print( wts.colptr(10*cols + 10));


        if ((yoff+blockysize) < rows)
        {
            yoff += linestowrite - rollback;
        }
        else
        {
            yoff=rows;
        }

        status = GDALTermProgress( yoff/(1.0*rows), NULL, NULL);
    }

    t_end = getWallTime();

    std::cout << "despeckamp processing time: " << (t_end-t_start)/60.0 << " mins \n";

    //Close the datasets
    GDALClose(inDataset);
    GDALClose(outDataset);
    GDALClose(wtsDataset);

    return(0);
       
};



/* Main driver*/
int main(int  argc, const char *argv[] ) {

    //Options  
    despeckOptions opts;
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
    status = despeck_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing calamp \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);
       
};
