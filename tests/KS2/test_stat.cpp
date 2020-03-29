#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include "KS2sample.hpp"


int main(int argc, char **argv)
{
   
    std::vector<std::string> files(4);
    files[0] = "Rayleigh1.bin";
    files[1] = "Rayleigh2.bin";
    files[2] = "Rayleigh3.bin";
    files[3] = "Exponential.bin";

    for(int ii=0; ii<4; ii++)
    {
        arma::fvec alldata;
        alldata.load(files[ii], arma::raw_binary);

        int nrows = alldata.size() / 2;
        arma::fmat x(alldata.memptr(), nrows,2, false);

        KS2sample worker(nrows);
        
        double prob = worker.test(x.colptr(0), x.colptr(1));
        std::cout << "File: " << files[ii] << " Prob: " << prob << "\n";
    }

}
