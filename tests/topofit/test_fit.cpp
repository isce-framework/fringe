#include <iostream>
#include <armadillo>
#include <vector>
#include <string>
#include "fringe_topofit.hpp"


int main(int argc, char **argv)
{
   
    std::vector<std::string> files(4);
    files[0] = "test1.dat";
    files[1] = "test2.dat";
    files[2] = "test3.dat";

    std::complex<float> cJ(0.0, 1.0);

    for(int ii=0; ii<3; ii++)
    {
        arma::fvec alldata;
        alldata.load(files[ii], arma::raw_binary);

        int nifg = (alldata.size() - 4)/2;
        arma::fvec bperp(alldata.memptr() + 4, nifg, false);
        arma::fvec ph(alldata.memptr()+nifg+4, nifg, false);
        
        arma::cx_fvec cph(nifg);
        arma::cx_fvec resid(nifg);

        for(int kk=0; kk<nifg; kk++)
            cph[kk] = std::exp(cJ * ph[kk]);

        float wvl = alldata[0];
        float rng = alldata[1];
        float inc = alldata[2];
        float delz = alldata[3];

        float Kmod, Cph;

        topofit worker(20.0, 0.1, wvl, nifg);
        
        float coh = worker.fit(cph.memptr(), bperp.memptr(), rng, inc, resid.memptr(), Kmod, Cph);
        std::cout << "File: " << files[ii] << " K: " << Kmod
                << " C: " << Cph << " coh: " << coh 
                << " K_true: " << delz/(wvl*rng*sin(inc*M_PI/180.0)/4.0/M_PI)<< "\n";
    }
}
