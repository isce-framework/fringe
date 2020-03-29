#ifndef FRINGE_DESPECK_H
#define FRINGE_DESPECK_H

/***
 * \author Piyush Agram.
 *  */

#include <iostream> 
#include <string>
#include <vector>  
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>  
#include "fringe/args.hxx"
#include "fringe/ulongmask.hpp"

//Input options for despeck
struct despeckOptions
{
    //Inputs
    std::string inputDS;     //Input VRT with SLCs as bands
    std::string wtsDS;       //Weights to use for averaging

    //Outputs
    std::string outputDS;       //Mask dataset 

    //Parameters
    int blocksize;           //Block size in rows
    int memsize;             //Memory in MB that can be used by the program

    int ibands[2];           //Bands to use for ifg

    int Nx;   //Half window size in pixels
    int Ny;   //Hald window size in lines

    bool computeCoherence; //Compute coherence for ifgs

    //Functions
    despeckOptions();        //Default constructor 
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report
};


//Default constructor
despeckOptions::despeckOptions()
{
    blocksize = 64;
    memsize = 512;
    Nx = 5;
    Ny = 5;
    ibands[0] = 1;   //Master band
    ibands[1] = -1;   //Slave band
    computeCoherence = false;
}


//Command line parser
int despeckOptions::initFromCmdLine(int argc, const char **argv)
{

    args::ArgumentParser parser("Despeckle SLC magnitude/ IFG using a neighborhood mask");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> input(parser, "inputDS*", "Input dataset", {'i'});
    args::ValueFlag<std::string> wts(parser, "wgtsDS*", "Input neighbor mask", {'m'});
    args::ValueFlag<std::string> output(parser, "outputDS*", "Despeckled amplitude dataset", {'o'});
    args::ValueFlag<int> lpb(parser,"linesperblock", "Lines per block for processing", {'l'});
    args::ValueFlag<int> mems(parser,"memorysize", "Memory in Mb", {'r'});
    args::ValueFlag<int> xx(parser, "xsize", "Half window in x", {'x'});
    args::ValueFlag<int> yy(parser, "ysize", "Half window in y", {'y'});
    args::ValueFlagList<int> bb(parser, "band", "Band index to use in SLC / IFG", {'b'});
    args::Flag coh(parser, "coherence", "Replace amplitude with coherence", {'c'});

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
    catch(args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }


    if ( (!input) || (!output) || (!wts))
    {
        std::cerr << parser;
        return 1;
    }

    inputDS = args::get(input);
    wtsDS = args::get(wts);
    outputDS = args::get(output);

    if (xx) { Nx = args::get(xx);}
    if (yy) {Ny = args::get(yy);}
    if (lpb) {blocksize = args::get(lpb);}
    if (mems) {memsize = args::get(mems);}
    if (coh) {computeCoherence = args::get(coh);}
    if (bb)
    { 
        std::vector<int> inbands = args::get(bb);
        if ((inbands.size() > 2) || (inbands.size() == 0))
        {
            std::cout << "Atleast 1 /2 band numbers must be provided for analysis. \n";
            std::cout << "Error parsing command line. Exiting with error code ... \n";
            return 1;
        }
        
        ibands[0] = inbands[0];
        if (inbands.size() == 2)
        {
            ibands[1] = inbands[1];
        }
        else
        {
            computeCoherence = false;
        }
        
    }

    return 0;
}


//Print report
void despeckOptions::print()
{

    std::cout << "Input Dataset: " << inputDS << std::endl;
    std::cout << "Weights Dataset: " << wtsDS << std::endl;
    std::cout << "Output Dataset: " << outputDS << std::endl;
    std::cout << "Window size: " << Nx << " " << Ny << std::endl;
    std::cout << "Memsize: " << memsize << " Mb \n";
    std::cout << "Blocksize: " << blocksize << " lines \n";
    std::cout << "Master band: " << ibands[0] << "\n";

    if (ibands[1] > 0)
    {
        std::cout << "Slave band: " << ibands[1] << "\n";
        if (computeCoherence)
            std::cout << "Coherence requested instead of amplitude \n";
    }
}
#endif //FRINGE_DESPECK_H
