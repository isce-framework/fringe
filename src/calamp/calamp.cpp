/**\file calamp.c
 * \author Piyush Agram.
 * \brief Amplitude calibration of SLCs.
 *  */

#include "calamp.hpp"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "fringe/fringe_common.hpp"

int calamp_process(calampOptions *opts)
{

    //Print user options to screen
    opts->print();


    int cols, rows, nbands;
    int blockysize;


    //First thing to do is to read dimenions.
    GDALDataset* inDataset = NULL;
    GDALDataset* mskDataset = NULL;
    GDALDataset* outDataset = NULL;

    //Clock variables
    double t_start, t_end;

    //Register GDAL drivers
    GDALAllRegister();


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

    //Get the dimensions and number of bands
    cols = inDataset->GetRasterXSize();
    rows = inDataset->GetRasterYSize();
    nbands = inDataset->GetRasterCount();

    std::cout << "Number of rows  = " << rows << "\n";
    std::cout << "Number of cols  = " << cols << "\n";
    std::cout << "Number of bands = " << nbands << "\n";


    //Check for mask dataset
    if (! opts->maskDS.empty())
    {
        if (strcasecmp(opts->maskDS.c_str(), "None") != 0)  
        {
            mskDataset = reinterpret_cast<GDALDataset *>( GDALOpenShared( opts->maskDS.c_str(), GA_ReadOnly));
            
            if (mskDataset == NULL)
            {
                std::cout << "Cannot open mask file { " << opts->maskDS << " } with GDAL for reading. \n";
                std::cout << "Exiting with error code .... (103) \n";
                GDALDestroyDriverManager();
                return 103;
            }

            //Read the dimensions of mask data to confirm its of same size
            int cc = mskDataset->GetRasterXSize();
            int rr = mskDataset->GetRasterYSize();
            int bb = mskDataset->GetRasterCount();
            int code = 0;

            if (cc != cols)
            {
                std::cout << "Mask file width does not match stack size width \n";
                code = 104;
            }
            if (rr != rows)
            {
                std::cout << "Mask file length does not match stack size length \n";
                code = 105;
            }
            if (bb != 1)
            {
                std::cout << "Mask file has more than one band \n";
                code = 106;
            }

            if (code !=0 )
            {
                std::cout << "Exiting with error code .... (" << code <<")\n";
                GDALClose(inDataset);
                GDALClose(mskDataset);
                GDALDestroyDriverManager();
                return code;
            }
        }
    }

    //Determine blocksizes
    blockysize = opts->blocksize * (int((opts->memsize * 1.0e6)/cols) / (opts->blocksize * 8));
    if (blockysize < opts->blocksize)
        blockysize = opts->blocksize;

    //If image is smaller than blocksize
    if (blockysize > rows)
        blockysize = rows;

    std::cout << "Block size = " << blockysize << " lines \n";

    //Start the clock
    t_start = getWallTime();

    //Array for storing SLC data that is read in
    std::vector< std::complex<float> > cpxdata(blockysize * cols);

    //Array for storing mask data that is read in 
    std::vector< unsigned char > maskdata( cols * blockysize, 1);

    //Array for storing sums for each band of data
    std::vector< double> totalsum(nbands, 0.0);

    //Array for storing normalization factors for each band of data
    std::vector< double> norms(nbands, 0.0);

    //Set up spacing requirements for data to be read in as BSQ 
    GSpacing pixelspacing = sizeof( std::complex<float>);
    GSpacing linespacing = pixelspacing * cols;
    GSpacing maskpixelspacing = sizeof(unsigned char);
    GSpacing masklinespacing = maskpixelspacing * cols;

    int totalblocks = ceil( rows / blockysize);
    int blockcount = 0;
    int status;
    int nthreads = 0;


    std::cout << "Total number of blocks to process: " << totalblocks << "\n";


    nthreads = numberOfThreads();

    //Block-by-block processing of data
    for (int yoff=0; yoff < rows; yoff += blockysize)
    {
        //Increment block counter
        blockcount++;

        //Determine number of lines to read
        int inysize = blockysize;
        if ((yoff+inysize) > rows)
            inysize = rows - yoff;


        //Read in a block of the mask
        if (mskDataset != NULL)
        {
            status = mskDataset->RasterIO(GF_Read, 0, yoff,
                                cols, inysize,
                                (void*) (&maskdata[0]),
                                cols, inysize, GDT_Byte,
                                1, NULL, maskpixelspacing,
                                masklinespacing, 0,  NULL);

            if (status != 0)
            {
                std::cout << "Error reading mask at line " << yoff << "\n";
                std::cout << "Exiting with error code .... (107) \n";
                GDALClose(inDataset);
                GDALClose(mskDataset);
                GDALDestroyDriverManager();
                return 107;
            }
        }
        
        //For each SLC
        for(int bb=0; bb < nbands; bb++)
        {
            //Read in the SLC data
            status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,   
                                cols, inysize,
                                (void*) (&cpxdata[0]), 
                                cols, inysize, GDT_CFloat32,
                                pixelspacing, linespacing, NULL);

            if (status != 0)
            {
                std::cout << "Error reading data from band " << bb+1 << " at line " 
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(inDataset);
                if (mskDataset != NULL)
                    GDALClose(mskDataset);
                GDALDestroyDriverManager();
                return 108;
            }



            double absval;
            int valid;
            double blocksum = 0.0;
            double blocknorm = 0.0;

            //Evaluate the sums of data and valid pixels
    #pragma omp parallel for \
            default(shared) \
            private(absval, valid) \
            reduction(+:blocksum,blocknorm) 
            for(int ii=0; ii< (cols * inysize); ii++)
            {
                absval = std::abs(cpxdata[ii]);
                absval = std::isnan(absval)? 0.0 : absval;
                valid = (absval != 0.0) * (maskdata[ii] > 0) ;
                blocksum += valid * absval;
                blocknorm += valid;
            }

            //Update the image sums with the blocksums
            totalsum[bb] += blocksum;
            norms[bb] += blocknorm;

            status = GDALTermProgress( ((blockcount-1)*nbands+bb)/(1.0*nbands*totalblocks), NULL, NULL);
        }
    }

    t_end = getWallTime();

    std::cout << "Calamp processing time: " << (t_end-t_start)/60.0 << " mins \n";

    std::cout << "Normalization coefficients \n";
    for(int ii=0; ii<nbands; ii++)
    {
        if (norms[ii] == 0)
        {
            std::cout << "No valid data found in Band " << ii+1 << ". Set to 1.0 \n";
            totalsum[ii] = 1.0;
        }
        else
        {
            totalsum[ii] /= norms[ii];
            std::cout << "Band " << ii+1 << ":" << totalsum[ii] << "  from  " <<  int(10000 * norms[ii]/ (1.0*rows*cols))/100.0 <<  " % of image \n";
        }
    }


    //Create output VRT
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("VRT");
    outDataset = (GDALDataset*) poDriver->CreateCopy(opts->outputDS.c_str(), inDataset, false,NULL, NULL, NULL); 
    
    for (int bb=1; bb<=nbands; bb++)
    {
        std::ostringstream strs;
        strs << totalsum[bb-1];
        std::string strval = strs.str();

        GDALRasterBand* bnd = outDataset->GetRasterBand(bb);
        GDALSetMetadataItem(bnd, "amplitudeConstant", strval.c_str(),  "slc");
    }


    //Close the datasets
    GDALClose(inDataset);
    GDALClose(outDataset);
    if (mskDataset != NULL)
        GDALClose(mskDataset);


  return(0);
       
};



/* Main driver*/
int main(int  argc, const char *argv[] ) {

    //Options  
    calampOptions opts;
    int status;

    //Parse command line options
    status = opts.initFromCmdLine(argc, argv);
    if (status != 0)  //Help message
    {
        std::cout << "Error processing command line for calamp \n";
        std::cout << "Exiting with non-zero return code .... \n";
        return(status);
    }

    //Execute 
    status = calamp_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing calamp \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);
       
};
