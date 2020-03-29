#ifndef FRINGE_KS2SAMPLE_H
#define FRINGE_KS2SAMPLE_H

#include <cmath>


/*Constants from ROOT. 
 * Source: ROOT/math/mathcore/src/Tmath.cxx: KolmogorovProb*/
struct KS2sample
{
    const double fj[4]={-2,-8,-18,-32};
    const double w = 2.50662827;

    // c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
    const double c1 = -1.2337005501361697;
    const double c2 = -11.103304951225528;
    const double c3 = -30.842513753404244;

    //Number of elements
    int length;

    KS2sample(){};
    KS2sample(int a){length = a;}
    double KolmogorovProb(const double z) const;
    double test(const float* a, const float* b) const;
};




/**********************************************************
 ****** Kolmogorov Smirnov 2-Sample test functions ********
 **********************************************************/


/*Kolmogorov Probability function
 * This function evaluates the Kolmogorov Probability function
 * "z" - Positive double 
 *
 * Source: ROOT/math/mathcore/src/Tmath.cxx: KolmogorovProb
 * All credit to ROOT developers / contributors */
inline double KS2sample::KolmogorovProb(const double z) const
{

    double u = std::abs(z);
    double p;
    double r[4];

    if (u < 0.2)
    {
        p = 1;
    }
    else if (u < 0.755)
    {
        double v = 1.0/(u*u);
        p = 1 - w*(::exp(c1*v) + ::exp(c2*v) + ::exp(c3*v))/u;
    }
    else if (u < 6.8116)
    {
        r[1] = 0;
        r[2] = 0;
        r[3] = 0;
        double v = u*u;
        int maxj = std::max(1, (int) (::round(3.0/u)));
        for(int j=0; j< maxj; j++)
        {
            r[j] = ::exp(fj[j]*v);
        }
        p = 2*(r[0] - r[1] + r[2] - r[3]);
    }
    else
    {
        p = 0;
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

inline double KS2sample::test(const float *a, const float *b) const
{
    int Na = length;
    int Nb = length;
    double rna = Na;
    double rnb = Nb;
    double sa = 1.0/rna;
    double sb = 1.0/rnb;

    int ia = 0;
    int ib = 0;
    double rdiff = 0;
    double rdmax = 0;

    for(int i=0; i< (Na+Nb); i++)
    {
        if(a[ia] < b[ib])
        {
            rdiff -= sa;
            ia++;
        }
        else if (a[ia] > b[ib])
        {
            rdiff += sb;
            ib++;
        }
        else
        {
            float x = a[ia];
            while(a[ia] == x && ia < Na)
            {
                rdiff -= sa;
                ia++;
            }
            while(b[ib] == x && ib < Nb)
            {
                rdiff += sb;
                ib++;
            }
        }

        if ((ia>=Na) || (ib>=Nb))
            break;

        rdmax = std::max(rdmax, std::abs(rdiff));
    }

    rdmax = std::max(rdmax, std::abs(rdiff));
    double z = rdmax * ::sqrt(rna*rnb/(rna+rnb));
//    std::cout << "K-stat: " << rdmax << "\n";
    double prob = KolmogorovProb(z);
    
    return prob;
}

#endif //FRINGE_KS2SAMPLE_H
