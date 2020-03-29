#ifndef FRINGE_CALAMP_H
#define FRINGE_CALAMP_H

/***
 * \author Piyush Agram.
 * \brief Amplitude calibration of interferograms.
 *
 *  Computes the amplitude calibration constant starting from interferograms directly.This is an ad-hoc technique directly related to computation
 *  of the amplitude calibration constants for coregistered SLCS like in the original StaMPS package. We use the square root of the absolute values
 *  of the amplitudes of the interferograms instead. Based on calamp.c written by Andy Hooper for StaMPS.The calibration constant for interferogram
 *   ''i'' is given by
 *
 * \f$ \mbox{calib}_i = (\displaystyle \sum \limits_{Valid pixels} |\mbox{Amp}|) \div \mbox{num pixels}\f$
 *  */

#include <iostream> 
#include <string.h>
#include <vector>  
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>     
#include "fringe/args.hxx"

//Input options for calamp
struct calampOptions
{
    //Inputs
    std::string inputDS;    //Input VRT with SLCs as bands
    std::string maskDS;     //Optional input GDAL dataset with mask

    //Outputs
    std::string outputDS;   //Output VRT with updated metadata including calamp coeffs

    //Parameters
    double defaultValue;    //Default valye to use for calamp coeff
    bool applySqrt;         //Apply sqrt when estimating constant
    int blocksize;          //Block size in rows
    int memsize;            //Memory in MB that can be used by the program

    //Functions
    calampOptions();        //Default constructor 
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report
};


//Default constructor
calampOptions::calampOptions()
{
    defaultValue = 1.0;
    applySqrt = false;
    blocksize = 128;
    memsize = 256;
}

int calampOptions::initFromCmdLine(int argc, const char **argv)
{
    args::ArgumentParser parser("Amplitude calibration using magnitudes of SLCS.");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> input(parser, "inputDS*", "Input Stack VRT", {'i'});
    args::ValueFlag<std::string> output(parser, "outputDS*", "Output Stack VRT", {'o'});
    args::ValueFlag<std::string> mask(parser, "maskDS", "Mask dataset", {'m'});
    args::ValueFlag<int> linesperblock(parser, "linesperblock", "Lines to load per block", {'l'});
    args::ValueFlag<int> memory(parser, "memoryinMb", "Memory to use in Mb", {'r'});
    args::ValueFlag<double> defval(parser, "defaultValue", "Default norm", {'n'}); 

    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 1;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    //Check mandatory inputs 
    if ((!input) || (!output))
    {
        std::cerr << parser;
        return 1;
    }
    
    inputDS = args::get(input);
    outputDS = args::get(output);

    if (mask) { maskDS = args::get(mask);}
    if (linesperblock) { blocksize = args::get(linesperblock);}
    if (memory) { memsize = args::get(memory);}
    if (defval) { defaultValue = args::get(defval);}

    return 0;
}


//Print report
void calampOptions::print()
{

    std::cout << "Input Dataset: " << inputDS << std::endl;
    std::cout << "Output Dataset: " << outputDS << std::endl;

    if (!maskDS.empty())
        std::cout << "Mask dataset: " << maskDS << std::endl;
    else
        std::cout << "No mask dataset has been defined." << std::endl;

    std::cout << "Lines per block: " << blocksize << std::endl;
    std::cout << "Memory size: " << memsize << " Mb \n";
    std::cout << "Default norm: " << defaultValue << "\n";

}
#endif //FRINGE_CALAMP_H
