#ifndef FRINGE_EVD_H
#define FRINGE_EVD_H

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

//Input options for evd
struct evdOptions
{
    //Inputs
    std::string inputDS;     //Input VRT with SLCs as bands
    std::string wtsDS;       //Self similar neighbot bit mask

    //Outputs
    std::string compSlc;     // name of compressed SLC
    std::string outputCompressedSlcFolder;  //output folder to store compressed SLCs
    std::string outputFolder;       //output stack
    std::string coherence;          //Coherence of the EVD description
    std::string compSLC;            //Compressed SLC

    int blocksize;           //Block size in rows
    int memsize;             //Memory in MB that can be used by the program
//    std::vector<int> ignorebands; //Bands to ignore for computing dispersion

    int Nx;   //Half window size in pixels
    int Ny;   //Hald window size in lines
    int minNeighbors; //Minimum number of neighbors

    int miniStackCount; //A counter for miniStack, used when Sequential algorithm is used. For one single stack of all acquisitions, miniStackCount = 1

    //Functions
    evdOptions();        //Default constructor 
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report
};


//Default constructor
evdOptions::evdOptions()
{
    blocksize = 64;
    memsize = 2048;
    Nx = 5;
    Ny = 5;
    minNeighbors = 2;
    miniStackCount = 1;
}


//Command line parser
int evdOptions::initFromCmdLine(int argc, const char **argv)
{

    args::ArgumentParser parser("Eigen value decomposition of a stack of SLCs");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> input(parser, "inputDS*", "Input Stack VRT", {'i'});
    args::ValueFlag<std::string> mask(parser, "wtsDS*", "Input neighbor mask", {'w'});
    args::ValueFlag<std::string> outfolder(parser, "output*", "Output folder for new stack", {'o'});
    args::ValueFlag<std::string> compfolder(parser, "outputComp*", "Output folder for compressed SLC", {'c'});
    args::ValueFlag<std::string> compname(parser, "compSlcName*", "Output compressed SLC", {'s'});
    args::ValueFlag<int> lpb(parser,"linesperblock", "Lines per block for processing", {'l'});
    args::ValueFlag<int> mems(parser,"memorysize", "Memory in Mb", {'r'});
    args::ValueFlag<int> xx(parser, "xsize", "Half window in x", {'x'});
    args::ValueFlag<int> yy(parser, "ysize", "Half window in y", {'y'});
    args::ValueFlag<int> mm(parser, "minNeighbors", "Minimum number of neighbors", {'n'});
    args::ValueFlag<int> kk(parser, "miniStackCount", "Mini-stack count", {'k'});
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


    if ( (!input) || (!mask) || (!outfolder))
    {
        std::cerr << parser;
        return 1;
    }

    inputDS = args::get(input);
    wtsDS = args::get(mask);
    outputFolder = args::get(outfolder);
    if (xx) { Nx = args::get(xx);}
    if (yy) {Ny = args::get(yy);}
    if (lpb) {blocksize = args::get(lpb);}
    if (mems) {memsize = args::get(mems);}
    if (mm) {minNeighbors = args::get(mm);}
    if (kk) {miniStackCount = args::get(kk);}
    if (compfolder) {
       outputCompressedSlcFolder = args::get(compfolder);
       if (compname) {
           compSLC = args::get(compname);
       }
       else{
           compSLC = "compslc.bin";
       }
    }
    else{
       outputCompressedSlcFolder = args::get(outfolder);
       compSLC = "compslc.bin";
    }

    return 0;
}


//Print report
void evdOptions::print()
{

    std::cout << "Input Dataset: " << inputDS << std::endl;
    std::cout << "Weights Dataset: " << wtsDS << std::endl;
    
    std::cout << "Output Folder " << outputFolder << std::endl;
    std::cout << "Output compressed SLC folder " << outputCompressedSlcFolder << std::endl;
    std::cout << "Output compressed SLC " << compSlc << std::endl;
    std::cout << "Window size: " << Nx << " " << Ny << std::endl;
    std::cout << "Memsize: " << memsize << " Mb \n";
    std::cout << "Blocksize: " << blocksize << " lines \n";
    std::cout << "Minimum neighbors: " << minNeighbors << "\n";
    std::cout << "Mini-stack counter: " << miniStackCount << "\n";
}   
#endif //FRINGE_EVD_H
