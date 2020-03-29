#ifndef FRINGE_TOPOFIT_H
#define FRINGE_TOPOFIT_H

#include <iostream>
#include <complex>
#include <algorithm>
#include <math.h>


struct topofit
{
    float maxError;    //Maximum DEM error
    float errorStep;   //Steps to search DEM error
    float wavelength;  //Radar wavelength
    int nHeights;        //Number of data points to fit.
    int nPoints;        //Number of interferograms for fit.
    int nSecHeights;     //Secondary search at 1/10 of errorStep;

    std::complex<float>  *cxworker;  //Temporary arrays
    float* fworker1;    


    //Functions
    topofit(const float, const float, const float, const int);
    ~topofit();
    inline float fit(const std::complex<float>*, const float *, 
            const float, const float, 
            std::complex<float>*, float& , float&);
};

//Constructor
topofit::topofit(const float maxErr , const float delta, const float wvl, const int nifg)
{
    maxError = maxErr;
    errorStep = delta;
    wavelength = wvl;
    nPoints = nifg;
    nHeights = (int) (2 * maxErr/delta + 1);
    nSecHeights = 21;

    cxworker = new std::complex<float>[nPoints];

    int narr = std::max(nHeights, nSecHeights);
    fworker1 = new float[narr];

}

topofit::~topofit()
{
    delete [] cxworker;
    delete [] fworker1;
}

inline float topofit::fit(const std::complex<float>* cpxphase, const float* bperp,
        const float rng, const float inc,
        std::complex<float> *resid, float& Kout, float& Cph)
{
    //Coherence will be returned at the end of the function
    float coherence;
    
    //Scale factor - constant for given geometry 
    float K = 4 * M_PI / (wavelength * rng * sin( M_PI * inc/180.0));

    //For each height candidate
    for(int ii=0; ii < nHeights; ii++)
    {
        std::complex<float> tempe(0,0);
        float tempf=0;
        float Kmod = K * (-maxError + ii * errorStep);

        for(int jj=0; jj<nPoints; jj++)
        {
            std::complex<float> tempc(0,Kmod*bperp[jj]);
            tempe += cpxphase[jj] * std::exp(tempc);
            tempf += abs(cpxphase[jj]);
        }
        fworker1[ii] = abs(tempe)/tempf;

//        std::cout << " ii = " << ii << " Kmod = " << Kmod << " coh = " << fworker1[ii] << "\n";
    }  

    //Preliminary height estimate
    int argmax = std::distance(fworker1, std::max_element(fworker1, fworker1 + nHeights));
    float estH = -maxError + argmax * errorStep;

//    std::cout << "Prelim K_est: " << K*estH << "  " << fworker1[argmax] << " arg: " << argmax << "\n";
//    std::cout << "Left: " << K*(estH-errorStep) << " " << fworker1[argmax-1] << "\n";
//    std::cout << "Right: " << K*(estH + errorStep) << " "<< fworker1[argmax+1] << "\n";
    
    
    //Start search around preliminary estimate
    float delhgt = 2*errorStep / (nSecHeights-1.0);
//    std::cout << "Steps : " << errorStep << " " << delhgt << "\n";
    for (int ii=0; ii < nSecHeights; ii++)
    {
        std::complex<float> tempe(0,0);
        float tempf = 0;
        float Kmod = K * (estH  - errorStep + ii * delhgt);

        for (int jj = 0; jj < nPoints; jj++)
        {
            std::complex<float> tempc(0, Kmod * bperp[jj]);
            tempe += cpxphase[jj] * std::exp(tempc);
            tempf += abs(cpxphase[jj]);
        }

        fworker1[ii] = abs(tempe)/tempf;
//        std::cout << "ii = " << ii <<  " Kmod = " << Kmod << " coh = " << fworker1[ii] << "\n";
    }
   
    //Final estimate of residual
    argmax = std::distance(fworker1, std::max_element(fworker1, fworker1 + nSecHeights));
    estH += (-errorStep + argmax * delhgt);

//    std::cout << "Final K_est: " << K*estH << " " << fworker1[argmax]  << " arg: " << argmax << "\n";

    {
        float Kmod = K*estH;
        std::complex<float> tempe(0,0);
        float tempf = 0;
        for(int jj=0; jj < nPoints; jj++)
        {
            std::complex<float> tempc(0, Kmod * bperp[jj]);
            resid[jj] = cpxphase[jj] * std::exp(tempc);
            tempe += resid[jj];
            tempf += abs(resid[jj]);
//              std::cout << jj << " " << std::abs(cpxphase[jj]) << " " << std::abs(resid[jj]) << "\n";
        }

        //Pack coherence into real part
        //Pack Kmod into imaginary part
        coherence = std::abs(tempe) / tempf;
        Cph = std::arg(tempe);
        Kout = K*estH;
    }


    return coherence;
}



#endif //FRINGE_TOPOFIT_H
