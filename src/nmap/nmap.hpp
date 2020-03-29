#ifndef FRINGE_NMAP_H
#define FRINGE_NMAP_H

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
#include "fringe/fringe_common.hpp"
#include "KS2sample.hpp"
#include "AD2unique.hpp"
#include "nmap_cuda.h"

//Input options for nmap
struct nmapOptions
{
    //Inputs
    std::string inputDS;     //Input VRT with SLCs as bands
    std::string maskDS;      //Mask dataset
    bool noGPU;             //Use GPU if available

    //Outputs
    std::string wtsDS;        //Neighborhood bit mask
    std::string ncountDS;    //Neighbor count

    //Parameters
    std::string method;       //KS2 or AD2

    int blocksize;           //Block size in rows
    int memsize;             //Memory in MB that can be used by the program
    std::vector<int> ignorebands; //Bands to ignore for computing dispersion

    int Nx;   //Half window size in pixels
    int Ny;   //Hald window size in lines

    double prob; //Minimum probability threshold 

    //Functions
    nmapOptions();        //Default constructor 
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report
};

//Function declaration
//int nmap_process(nmapOptions *opt);

//Default constructor
nmapOptions::nmapOptions()
{
    blocksize = 64;
    memsize = 512;
    Nx = 5;
    Ny = 5;
    method = "KS2";
    prob = 0.05;
    noGPU = false;
}


//Command line parser
int nmapOptions::initFromCmdLine(int argc, const char **argv)
{

    args::ArgumentParser parser("Create neighborhood mask and count map using KS statistics");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> input(parser, "inputDS*", "Input Stack VRT", {'i'});
    args::ValueFlag<std::string> output(parser, "outputDS*", "Output neighbor mask", {'o'});
    args::ValueFlag<std::string> count(parser, "countDS*", "Neighbor count dataset", {'c'});
    args::ValueFlag<std::string> mask(parser, "maskDS*", "Mask dataset", {'m'});
    args::ValueFlag<int> lpb(parser,"linesperblock", "Lines per block for processing", {'l'});
    args::ValueFlag<int> mems(parser,"memorysize", "Memory in Mb", {'r'});
    args::ValueFlag<int> xx(parser, "xsize", "Half window in x", {'x'});
    args::ValueFlag<int> yy(parser, "ysize", "Half window in y", {'y'});
    args::ValueFlag<double> P(parser, "prob", "P-value threshold", {'p'});
    args::ValueFlag<std::string> tmethod(parser, "testMethod", "Statistics test - KS2 or AD2", {'t'});
    args::Flag nocuda(parser, "nogpu", "Do not use GPU", {"nogpu"});


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


    if ( (!input) || (!output) || (!count))
    {
        std::cerr << parser;
        return 1;
    }

    inputDS = args::get(input);
    wtsDS = args::get(output);
    ncountDS = args::get(count);

    if (mask) { maskDS = args::get(mask);}
    if (xx) { Nx = args::get(xx);}
    if (yy) {Ny = args::get(yy);}
    if (P)  {prob = args::get(P);}
    if (lpb) {blocksize = args::get(lpb);}
    if (mems) {memsize = args::get(mems);} 
    if (tmethod)
    {
        std::string inmethod = args::get(tmethod);
        if ( (inmethod.compare("KS2")==0) || (inmethod.compare("AD2")==0))
        {
            method = inmethod;
        }
        else
        {
            std::cout << "Input method should either be KS2 or AD2 \n";
            std::cerr << parser;
            return 1;
        }
    }

    if (nocuda)
    {
        noGPU = true;
    }
    return 0;
}


//Print report
void nmapOptions::print()
{

    std::cout << "Input Dataset: " << inputDS << std::endl;
    std::cout << "Weights Dataset: " << wtsDS << std::endl;
    std::cout << "Count Dataset: " << ncountDS << std::endl;
    std::cout << "Mask Dataset: " << maskDS << std::endl;
    std::cout << "Window size: " << Nx << " " << Ny << std::endl;
    std::cout << "Threshold : " << prob << std::endl;
    std::cout << "Memsize: " << memsize << " Mb \n";
    std::cout << "Blocksize: " << blocksize << " lines \n";
    std::cout << "GPU user request: " << !(noGPU) << " \n";

}
#endif //FRINGE_NMAP_H
