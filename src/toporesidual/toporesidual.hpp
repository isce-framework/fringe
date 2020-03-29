#ifndef FRINGE_TOPO_H
#define FRINGE_TOPO_H

/***
 *  \author Heresh Fattahi.
 *   */

/*#include <cmath>
#include <cstdio>
#include <cstdlib>
*/

#include <iostream>
#include <string.h>
#include <armadillo>
#include <vector>
#include <string>
#include <complex>
#include "gdal.h"
#include "gdal_priv.h"
#include "fringe_topofit.hpp"
#include "fringe/args.hxx"
#include "fringe/fringe_common.hpp"
struct topoOptions
{

    //Inputs
    std::string phaseDS;   //stack of cimplex interferograms or complex phase time-series
    std::string bperpDS;   // stack of perpendicular baseline 
    std::string incDS;     // incidence angle dataset
    
    //outputs
    std::string outputFolder;
    std::string cohDS="tempCoh.bin";     // temporal coherence dataset
    std::string deltazDS="demError.bin";   // estimated topographic residual

    int blocksize;           //Block size in rows
    int memsize;             //Memory in MB that can be used by the program
    
    float wavelength;

    topoOptions();        //Default constructor
    int initFromCmdLine(int, const char**); //Command line parser
    void print();   //Print Report        
    
};

topoOptions::topoOptions()
{
    blocksize = 64;
    memsize = 2048;
    cohDS = "toporesidual_tempCoh.coh";
    deltazDS = "toporesidual.dz";

}

int topoOptions::initFromCmdLine(int argc, const char **argv)
{
    args::ArgumentParser parser("Topographical residual estimation.");
    args::HelpFlag flag(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> phase(parser, "phaseDS*", "Input Stack VRT of the wrapped phase time-series", {'p'});
    args::ValueFlag<std::string> bperp(parser, "bperpDS*", "Input Stack VRT of the bperp time-series", {'b'});
    args::ValueFlag<std::string> inc(parser, "incDS*", "Input incidence angle dataset", {'i'});
    args::ValueFlag<std::string> output(parser, "outputFolder", "Output folder", {'o'});
    args::ValueFlag<std::string> coh(parser, "cohDS*", "Output temporal coherence", {'c'});
    args::ValueFlag<std::string> dz(parser, "dzDS*", "Output topographical residual dataset coherence", {'z'});
    args::ValueFlag<int> linesperblock(parser, "linesperblock", "Lines to load per block", {'l'});
    args::ValueFlag<int> memory(parser, "memoryinMb", "Memory to use in Mb", {'r'});

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

    if ((!phase) || (!bperp) || (!inc) || (!output))
    {
        std::cerr << parser;
        return 1;
    }

    phaseDS = args::get(phase);
    bperpDS = args::get(bperp);
    incDS = args::get(inc);
    outputFolder = args::get(output);

    if (coh) {cohDS = args::get(coh);}
    if (dz) {deltazDS = args::get(dz);}
    if (linesperblock) { blocksize = args::get(linesperblock);}
    if (memory) { memsize = args::get(memory);}

    return 0;

}

//Print report
void topoOptions::print()
{

    std::cout << "Input wrapped phase Dataset: " << phaseDS << std::endl;
    std::cout << "Input baseline Dataset: " << bperpDS << std::endl;
    std::cout << "Input incidence angle Dataset: " << incDS << std::endl;
    std::cout << "Output estimated topographical residual dataset: " << deltazDS << std::endl;
    std::cout << "Output temporal coherence dataset: " << cohDS << std::endl;
    std::cout << "Lines per block: " << blocksize << std::endl;
    std::cout << "Memory size: " << memsize << " Mb \n";

}


#endif




