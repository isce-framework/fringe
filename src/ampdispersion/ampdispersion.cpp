/**\file ampdispersion.c
 * \author Piyush Agram.
 *  */

#include "ampdispersion.hpp"
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "fringe/fringe_common.hpp"

int ampdispersion_process(ampdispersionOptions *opts)
{

    //Print user options to screen
    opts->print();

    int cols, rows, nbands;
    int blockysize;


    //First thing to do is to read dimenions.
    GDALDataset* inDataset = NULL;
    GDALDataset* daDataset = NULL;
    GDALDataset* meanampDataset = NULL;

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

    //Determine the number of bands and dimensions
    cols = inDataset->GetRasterXSize();
    rows = inDataset->GetRasterYSize();
    nbands = inDataset->GetRasterCount();

    std::cout << "Number of rows  = " << rows << "\n";
    std::cout << "Number of cols  = " << cols << "\n";
    std::cout << "Number of bands = " << nbands << "\n";


    //Determine blocksizes
    blockysize = opts->blocksize * (int((opts->memsize * 1.0e6)/cols) / (opts->blocksize * 8) * 4);
    if (blockysize < opts->blocksize)
        blockysize = opts->blocksize;

    //If image is smaller than block size
    if (blockysize > rows)
        blockysize = rows;

    std::cout << "Block size = " << blockysize << " lines \n";

    int totalblocks = ceil( rows / (1.0 * blockysize));
    std::cout << "Total number of blocks to process: " << totalblocks << "\n";

    //Start the clock
    t_start = getWallTime();

    //Array to read in complex SLC data
    std::vector< std::complex<float> > cpxdata(blockysize * cols);

    //Array to store means for the pixels
    std::vector<double> mean(blockysize*cols);

    //Array to store means of squares of the pixels
    std::vector<double> meansq(blockysize*cols);

    //Array to store normalization values
    std::vector<double> norms(blockysize*cols);

    int blockcount = 0;
    int status;
    int nthreads = 0;
    std::vector<double> alpha(nbands);

    //Populate norms by reading in specific metadata 
    for (int bb=0; bb < nbands; bb++)
    {
        const char *bandconst = inDataset->GetRasterBand(bb+1)->GetMetadataItem( "amplitudeConstant", "slc");
        double bconst = 0.0;

        if (bandconst != NULL)
        {
             bconst = CPLScanDouble(bandconst, 20);
        }
        
        if(bconst <= 0.0)
        {
            bconst = 1.0;
            std::cout << "No calibration constant found for band " << bb+1 <<". Setting to 1.0. \n";
        }
        alpha[bb] = bconst;
    }


    if ((opts->refband < 1) ||(opts->refband > nbands))
    {
        std::cout <<"Reference band number: " << opts->refband << " is invalid \n";
        std::cout << "Exiting with non-zero error code .... (102) \n";
        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 102;
    }


    //Normalize w.r.t to a single reference band
    double normval = alpha[opts->refband-1];
    for (int bb=0; bb < nbands; bb++)
    {
            alpha[bb] /= normval;
    }
    alpha[opts->refband-1] = 1.0;

    for(int bb=0; bb < nbands; bb++)
        std::cout << "Band " << bb+1 << ": " << alpha[bb] << "\n";

   
    //Creating output dataset for storing amplitude dispersion 
    GDALDriver *poDriver = (GDALDriver*) GDALGetDriverByName("ENVI");
    char **mOptions = NULL;
    mOptions = CSLSetNameValue(mOptions, "INTERLEAVE", "BIP");
    mOptions = CSLSetNameValue(mOptions, "SUFFIX", "ADD");

    daDataset = (GDALDataset*) poDriver->Create(opts->daDS.c_str(), cols, rows, 1, GDT_Float32, mOptions);

    if (daDataset == NULL)
    {
        std::cout << "Could not create ampdisp dataset {" << opts->daDS << "} \n";
        std::cout << "Exiting with non-zero error code ... 103 \n";

        GDALClose(inDataset);
        GDALDestroyDriverManager();
        return 103;
    }

    //Creating output dataset for storing mean amplitude
    meanampDataset = (GDALDataset*) poDriver->Create(opts->meanampDS.c_str(), cols, rows, 1, GDT_Float32, mOptions); 
    if (meanampDataset == NULL)
    {
        std::cout << "Could not create meanamp dataset {" << opts->meanampDS << "} \n";
        std::cout << "Exiting with non-zero error code ... 104 \n";

        GDALClose(inDataset);
        GDALClose(daDataset);
        GDALDestroyDriverManager();
        return 104;
    }

    CSLDestroy(mOptions);

    nthreads = numberOfThreads();

    //Block-by-block processing
    for (int yoff=0; yoff < rows; yoff += blockysize)
    {
        //Increment block counter
        blockcount++;

        //Determine number of rows to read
        int inysize = blockysize;
        if ((yoff+inysize) > rows)
            inysize = rows - yoff;

        //Initialize sums to zero
        memset( &mean[0], 0, mean.size() * sizeof(mean[0]));
        memset( &meansq[0], 0, meansq.size() * sizeof(meansq[0]));
        memset( &norms[0], 0, norms.size() * sizeof(norms[0]));

        //For each SLC
        for(int bb=0; bb < nbands; bb++)
        {
            //Read in the SLC data
            status = inDataset->GetRasterBand(bb+1)->RasterIO( GF_Read, 0, yoff,   
                                cols, inysize,
                                (void*) (&cpxdata[0]), 
                                cols, inysize, GDT_CFloat32,
                                sizeof(std::complex<float>),
                                sizeof(std::complex<float>)*cols, NULL);

            if (status != 0)
            {
                std::cout << "Error reading data from band " << bb+1 << " at line " 
                    << yoff << "\n";
                std::cout << "Exiting with error code .... (108) \n";
                GDALClose(inDataset);
                GDALClose(daDataset);
                GDALClose(meanampDataset);
                GDALDestroyDriverManager();
                return 108;
            }
            
            //Update the sums with data from each SLC
    #pragma omp parallel for \
            default(shared) 
            for(int ii=0; ii< (cols * inysize); ii++)
            {
                double absval = std::abs(cpxdata[ii]);
                int valid = (absval != 0.0);
                absval *= (valid/alpha[bb]);
                mean[ii] += absval;
                meansq[ii] += absval*absval;
                norms[ii] += valid;
            }
        }

        //Perform the computation for ampdisperion
    #pragma omp parallel for \
        default(shared) 
        for(int ii=0; ii < (cols * inysize); ii++)
        {
            double avg, avg2, sdev;
            if (norms[ii] > 1)
            {
                avg = mean[ii]/norms[ii];
                avg2 = meansq[ii]/norms[ii];
                sdev = ::sqrt(avg2-avg*avg);
                mean[ii] = avg;
                meansq[ii] = ((!std::isnan(sdev)) && (sdev > 0))? (sdev/ avg) : -1;
            }
            else
            {
                mean[ii] = 0.0;
                meansq[ii] = -1.0;
            }
        }


        //Write result to amplitude dispersion dataset
        status = daDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff,
                    cols, inysize,
                    (void*) (&meansq[0]), 
                    cols, inysize, GDT_Float64,
                    sizeof(double), sizeof(double)*cols, NULL);

        if (status != 0)
        {
            std::cout << "Error write ampdisp data at line " << yoff << "\n";
            std::cout << "Exiting with error code .... (109) \n";
            GDALClose(inDataset);
            GDALClose(daDataset);
            GDALClose(meanampDataset);
            GDALDestroyDriverManager();
            return 109;
        }


        //Write result to mean amplitude dataset
        status = meanampDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, yoff,
                    cols, inysize,
                    (void*) (&mean[0]), 
                    cols, inysize, GDT_Float64,
                    sizeof(double), sizeof(double)*cols, NULL);

        if (status != 0)
        {
            std::cout << "Error write meanamp data at line " << yoff << "\n";
            std::cout << "Exiting with error code .... (110) \n";
            GDALClose(inDataset);
            GDALClose(daDataset);
            GDALClose(meanampDataset);
            GDALDestroyDriverManager();
            return 110;
        }

        status = GDALTermProgress( (yoff+inysize)/(1.0*rows), NULL, NULL);
    }

    t_end = getWallTime();

    std::cout << "ampdispersion processing time: " << (t_end-t_start)/60.0 << " mins \n";

    //Annotate the output datasets with number of SLCs/nbands used for future reference
    status = daDataset->SetMetadataItem("N", std::to_string(nbands).c_str(), "ENVI");
    status = meanampDataset->SetMetadataItem("N", std::to_string(nbands).c_str(), "ENVI");

    //Close the datasets
    GDALClose(inDataset);
    GDALClose(daDataset);
    GDALClose(meanampDataset);

    return(0);
       
};



/* Main driver*/
int main(int  argc, const char *argv[] ) {

    //Options  
    ampdispersionOptions opts;
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
    status = ampdispersion_process(&opts);

    if (status != 0)
    {
        std::cout << "Error processing calamp \n";
        std::cout << "Exiting with non-zero return code ....\n";
        return (status);
    }

    return(0);
       
};
