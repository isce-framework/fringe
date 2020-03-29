#ifndef FRINGE_AMPDISPERSION_H
#define FRINGE_AMPDISPERSION_H

/***
 * \author Piyush Agram.
 *  */

#include <iostream> 
#include <string.h>
#include <vector>  
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>     
#include "fringe/args.hxx"

//Input options for ampdispersion
struct ampdispersionOptions
{
    //Inputs
    std::string inputDS;     //Input VRT with SLCs as bands

    //Outputs
    std::string meanampDS;   //Mean normalized amplitude
    std::string daDS;        //Amplitude dispersion

    //Parameters
    int blocksize;           //Block size in rows
    int memsize;             //Memory in MB that can be used by the program
    std::vector<int> ignorebands; //Bands to ignore for computing dispersion
    int refband;             //Reference band - corresponding to master image

    //Functions
    ampdispersionOptions();        //Default constructor 
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report
};


//Default constructor
ampdispersionOptions::ampdispersionOptions()
{
    blocksize = 64;
    memsize = 256;
    refband = 1;
}


//Command line parser
int ampdispersionOptions::initFromCmdLine(int argc, const char **argv)
{
   
    args::ArgumentParser parser("Amplitude dispersion and mean amplitude calculator");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> input(parser, "inputDS*", "Input Stack VRT", {'i'});
    args::ValueFlag<std::string> output(parser, "dispDS*", "Output dispersion dataset", {'o'});
    args::ValueFlag<std::string> mean(parser, "meanDS*", "Mean amplitude dataset", {'m'});
    args::ValueFlag<int> lpb(parser, "linesperblock", "Lines to read per block", {'l'});
    args::ValueFlag<int> mems(parser, "memorysize", "Memory in Mb", {'s'});
    args::ValueFlag<int> ref(parser, "refband", "Reference band for norm", {'b'});

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

    if ((!input) || (!output) || (!mean))
    {
        std::cerr << parser;
        return 1;
    }

    inputDS = args::get(input);
    daDS = args::get(output);
    meanampDS = args::get(mean);


    if (lpb) { blocksize = args::get(lpb);}
    if (mems) { memsize = args::get(mems);}
    if (ref) { refband = args::get(ref);}


    return 0;
}




//Print report
void ampdispersionOptions::print()
{

    std::cout << "Input Dataset: " << inputDS << std::endl;
    std::cout << "Ampdisp Dataset: " << daDS << std::endl;
    std::cout << "Meanamp Dataset: " << meanampDS << std::endl;
    std::cout << "Lines per block: " << blocksize << std::endl;
    std::cout << "Memory size: " << memsize << std::endl;
    std::cout << "Ref band: " << refband << std::endl;

}
#endif //FRINGE_AMPDISPERSION_H
