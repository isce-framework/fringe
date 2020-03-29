#ifndef FRINGE_KS2SAMPLE_CUDA_H
#define FRINGE_KS2SAMPLE_CUDA_H

#include <cmath>
#include <cuda_runtime.h>

/*Constants from ROOT. 
 * Source: ROOT/math/mathcore/src/Tmath.cxx: KolmogorovProb*/

/**********************************************************
 ****** Kolmogorov Smirnov 2-Sample test functions ********
 **********************************************************/
/*Kolmogorov Probability function
 * This function evaluates the Kolmogorov Probability function
 * "z" - Positive double 
 *
 * Source: ROOT/math/mathcore/src/Tmath.cxx: KolmogorovProb
 * All credit to ROOT developers / contributors */
__device__ inline double KolmogorovProb(const double z) 
{

    double u = fabs(z);
    double p = 0.0;
    double r[4];
    int maxj;
    double v;
    if (u < 0.2)
    {
        p = 1;
    }
    else if (u < 0.755)
    {
        double v = 1.0/(u*u);
        p = 1 - 2.50662827*(exp(-1.2337005501361697*v) + exp(-11.103304951225528*v) + exp(-30.842513753404244*v))/u;
    }
    else if (u < 6.8116)
    {
        r[1] = 0;
        r[2] = 0;
        r[3] = 0;
        v = u*u;
        maxj = max(1, (int) (llround(3.0/u)));
        r[0] = exp(-2 * v);
        r[1] = (maxj > 1) * exp(-8*v);
        r[2] = (maxj > 2) * exp(-18*v);
        r[3] = (maxj > 3) * exp(-32*v);
        p = 2*(r[0] - r[1] + r[2] - r[3]);
    }

    return p;
}


/*Kolmogorov-Smirnov 2-sample test
 * Inputs are 2 arrays
 * array "a" of size "Na" 
 * array "b" of size "Nb" 
 *
 * Implementation adapted from CERN's ROOT Data Analysis Framework
 * https://root.cern.ch
 * All credit for this module goes to developers  contributors of the ROOT package.
 * Source: ROOT/math/mathcore/src/TMath.cxx:KolmogorovTest
 */

__device__ inline double KS2test(const float *a, const float *b, int length)
{
    int ia = 0;
    int ib = 0;
    double rdiff = 0;
    double rdmax = 0;
    double sa = 1.0/ (1.0 * length);
    float x;
    for(int i=0; i< (2*length); i++)
    {
        if(a[ia] < b[ib])
        {
            rdiff -= sa;
            ia++;
        }
        else if (a[ia] > b[ib])
        {
            rdiff += sa;
            ib++;
        }
        else
        {
            x = a[ia];
            while(a[ia] == x && ia < length)
            {
                rdiff -= sa;
                ia++;
            }
            while(b[ib] == x && ib < length)
            {
                rdiff += sa;
                ib++;
            }
        }

        if ((ia>=length) || (ib>=length))
            break;

        rdmax = fmax(rdmax, fabs(rdiff));
    }

    rdmax = fmax(rdmax, fabs(rdiff));
    sa = rdmax * sqrt(0.5 * length);
//    std::cout << "K-stat: " << rdmax << "\n";
    return KolmogorovProb(sa);
    
}

#endif //FRINGE_KS2SAMPLE_CUDA_H
