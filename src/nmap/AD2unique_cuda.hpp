#ifndef FRINGE_AD2UNIQUE_CUDA_H
#define FRINGE_AD2UNIQUE_CUDA_H

#include <cmath>

struct AD2unique
{

    __constant__ double ts[280] ={ -1.1954, -1.5806, -1.8172, 
                    -2.0032, -2.2526, -2.4204, -2.5283, -4.2649, -1.1786, -1.5394, 
                    -1.7728, -1.9426, -2.1685, -2.3288, -2.4374, -3.8906, -1.166, 
                    -1.5193, -1.7462, -1.9067, -2.126, -2.2818, -2.3926, -3.719, 
                    -1.1407, -1.4659, -1.671, -1.8105, -2.0048, -2.1356, -2.2348, 
                    -3.2905, -1.1253, -1.4371, -1.6314, -1.7619, -1.9396, -2.0637, 
                    -2.1521, -3.0902, -1.0777, -1.3503, -1.5102, -1.6177, -1.761, 
                    -1.8537, -1.9178, -2.5758, -1.0489, -1.2984, -1.4415, -1.5355, 
                    -1.6625, -1.738, -1.7936, -2.3263, -0.9978, -1.2098, -1.3251, 
                    -1.4007, -1.4977, -1.5555, -1.5941, -1.96, -0.9417, -1.1187, 
                    -1.209, -1.2671, -1.3382, -1.379, -1.405, -1.6449, -0.8981, -1.0491, 
                    -1.1235, -1.1692, -1.2249, -1.2552, -1.2755, -1.4395, -0.8598, 
                    -0.9904, -1.0513, -1.0879, -1.1317, -1.155, -1.1694, -1.2816, 
                    -0.7258, -0.7938, -0.8188, -0.8312, -0.8435, -0.8471, -0.8496, 
                    -0.8416, -0.5966, -0.617, -0.6177, -0.6139, -0.6073, -0.5987, 
                    -0.5941, -0.5244, -0.4572, -0.4383, -0.419, -0.4033, -0.3834, 
                    -0.3676, -0.3587, -0.2533, -0.2966, -0.2428, -0.2078, -0.1844, 
                    -0.1548, -0.1346, -0.1224, 0, -0.1009, -0.0169, 0.0304, 0.0596, 
                    0.0933, 0.1156, 0.1294, 0.2533, 0.1571, 0.2635, 0.3169, 0.348, 
                    0.3823, 0.4038, 0.4166, 0.5244, 0.5357, 0.6496, 0.6992, 0.7246, 
                    0.7528, 0.7683, 0.7771, 0.8416, 1.2255, 1.2989, 1.3202, 1.3254, 
                    1.3305, 1.3286, 1.3257, 1.2816, 1.5262, 1.5677, 1.5709, 1.5663, 
                    1.5561, 1.5449, 1.5356, 1.4395, 1.9633, 1.943, 1.919, 1.8975, 
                    1.8641, 1.8389, 1.8212, 1.6449, 2.7314, 2.5899, 2.5, 2.4451, 
                    2.3664, 2.3155, 2.2823, 1.96, 3.7825, 3.4425, 3.2582, 3.1423, 
                    3.0036, 2.9101, 2.8579, 2.3263, 4.1241, 3.716, 3.4984, 3.3651, 
                    3.2003, 3.0928, 3.0311, 2.4324, 4.6044, 4.0847, 3.8348, 3.6714, 
                    3.4721, 3.3453, 3.2777, 2.5758, 5.409, 4.7223, 4.4022, 4.1791, 
                    3.9357, 3.7809, 3.6963, 2.807, 6.4954, 5.5823, 5.1456, 4.8657, 
                    4.5506, 4.3275, 4.2228, 3.0902, 6.8279, 5.8282, 5.3658, 5.0749, 
                    4.7318, 4.4923, 4.3642, 3.1747, 7.2755, 6.197, 5.6715, 5.3642, 
                    4.9991, 4.7135, 4.5945, 3.2905, 8.1885, 6.8537, 6.2077, 5.8499, 
                    5.4246, 5.1137, 4.9555, 3.4808, 9.3061, 7.6592, 6.85, 6.4806, 
                    5.9919, 5.6122, 5.5136, 3.719, 9.6132, 7.9234, 7.1025, 6.6731, 
                    6.1549, 5.8217, 5.7345, 3.7911, 10.0989, 8.2395, 7.4326, 6.9567, 
                    6.3908, 6.011, 5.9566, 3.8906, 10.8825, 8.8994, 7.8934, 7.4501, 
                    6.9009, 6.4538, 6.2705, 4.0556, 11.8537, 9.5482, 8.5568, 8.0283, 
                    7.4418, 6.9524, 6.6195, 4.2649 };

      // p values bins 
    __constant__ double prob[35] = { .00001,.00005,.0001,.0005,.001,.005,.01,.025,.05,.075,.1,.2,.3,.4,.5,.6,.7,.8,.9,
                        .925,.95,.975,.99,.9925,.995,.9975,.999,.99925,.9995,.99975,.9999,.999925,
                        .99995,.999975,.99999 };

    __constant__ int ns = 8;
    __constant__ int nbins = 35;

    //Members
    int length;
    double sigmaNorm;
    double *ts2;
    double *lp;
    float *fCombinedSamples;
    float *z;
    int *fij;
    int *lvec;
    double *invI;


    AD2unique(){};
    AD2unique(int);
    ~AD2unique();

    void init(int);
    void prepare();
    double PValueADKSamples(const double) const;
    double getSigmaN(const int, const int);
    double test(const float*, const float*);

};


AD2unique::AD2unique(int num)
{
    init(num);
}

void AD2unique::init(int num)
{
    length = num;
    ts2 = new double[nbins];
    lp = new double[nbins];
    fCombinedSamples = new float[2*length];
    z = new float[2*length];
    fij = new int[4*length];
    lvec = new int[2*length];
    invI = new double[2*length];
    prepare();
}

AD2unique::~AD2unique()
{
    delete [] ts2;
    delete [] lp;
    delete [] fCombinedSamples;
    delete [] z;
    delete [] fij;
    delete [] lvec;
    delete [] invI;
}


inline void AD2unique::prepare()
{
    for(int i=0; i < nbins; i++)
    {
        ts2[i] = ts[i*ns];
        lp[i] = std::log( (1.0 - prob[i])/prob[i]);
    }

    sigmaNorm = getSigmaN(length, length);
}

/**********************************************************
 ******** Anderson Darling 2-Sample test functions ********
 **********************************************************/
inline double AD2unique::PValueADKSamples(const double tx) const
{

    int i1, i2;

    for (int ii=0; ii< nbins; ii++)
    {
        i1 = ii-1;
        if (tx <= ts2[ii])
            break;
    }
    i2 = i1+1;

    // if tx is before min of tabluated data
    if (i1 < 0) { 
         i1 = 0;
         i2 = 1;
    }
    
    // if tx is after max of tabulated data
    if (i2 >= nbins ) { 
         i1 = nbins-2; 
         i2 = nbins-1;
    }

    double lp1 = lp[i1]; 
    double lp2 = lp[i2];
    double tx1 = ts2[i1];
    double tx2 = ts2[i2];

    /// find interpolated (or extrapolated value)( 
    double lp0 = (lp1-lp2) * (tx - tx2)/ ( tx1-tx2) + lp2; 

    double p0 = exp(lp0)/(1. + exp(lp0) );
    return p0; 
}

inline double AD2unique::getSigmaN(const int ns1, const int ns2)
{
    double sigmaN = 0.0;
    double h = 0.0;
    double H = 0.0;
    double g = 0.0;
    double k=2.0;
    double a,b,c,d;

    int N = ns1 + ns2;
    H = 1.0/(1.0*ns1) + 1.0/(1.0*ns2);

    // use approximate formulas for large N
    // cache Sum( 1 / i)
    if (N < 2000)
    {
        for (int i = 1; i < N; ++i)
        {
            invI[i] = 1.0 / i;
            h += invI[i];
        }
        for (int i = 1; i < N - 1; ++i)
        {
            double tmp = invI[N-i];
            for (int j = i + 1; j < N; ++j)
            {
               g += tmp * invI[j];
            }
        }
    }
    else
    {
        const double emc = 0.5772156649015328606065120900824024; // Euler-Mascheroni constant
        h = std::log(double(N-1) ) + emc;
        g = (M_PI)*(M_PI)/6.0;
    }

    double k2 = std::pow(k,2);
    a = (4 * g - 6) * (k - 1)  + (10 - 6 * g) * H;
    b = (2 * g - 4) * k2 + 8 * h * k + (2 * g - 14 * h - 4) * H - 8 * h + 4 * g - 6;
    c = (6 * h + 2 * g - 2) * k2 + (4 * h - 4 *g + 6) * k + (2 * h - 6) * H + 4 * h;
    d = (2 * h + 6) * k2 - 4 * h * k;
    sigmaN +=  a * std::pow(double(N),3) + b * std::pow(double(N),2) + c * N + d;
    sigmaN /= ( double(N - 1) * double(N - 2) * double(N - 3) );
    sigmaN = sqrt(sigmaN);
    return sigmaN;
}


inline double AD2unique::test(const float *a, const float *b)
{

    int Na = length;
    int Nb = length;

    int ns1 = Na;
    int ns2 = Nb;
    int nsum = ns1 + ns2;


    {
        int ii=0;
        int jj=0;
        int kk=0;
        while( ii < ns1 && jj < ns2)
        {
            if (a[ii] < b[jj])
            {
                fCombinedSamples[kk] = a[ii];
                ii++;
                fij[2*kk] = 1;
                fij[2*kk+1] = 0;
            }
            else
            {
                fCombinedSamples[kk] = b[jj];
                jj++;
                fij[2*kk] = 0;
                fij[2*kk+1] = 1;
            }
            kk++;
        }

        if (ii>=ns1)
        {
            while(jj < ns2)
            {
                fCombinedSamples[kk] = b[jj];
                jj++;
                fij[2*kk] = 0;
                fij[2*kk+1] = 1;
                kk++;
            }
        }
        if (jj>=ns2)
        {
            while(ii < ns1)
            {
                fCombinedSamples[kk] = a[ii];
                ii++;
                fij[2*kk] = 1;
                fij[2*kk+1] = 0;
                kk++;
            }
        }
    }
    
    int l=2*length;

    //adkTest

    double mij;
    double maij;
    double innerSum;
    double aInnerSum;
    double bj;
    double baj;
    double tmp;

    double AkN2 = 0.0;
    double AakN2 = 0.0;

    //Start with array a
    mij = 0;
    maij = 0;
    innerSum = 0;
    aInnerSum = 0;

    for (int j=0; j<l; j++)
    {
        mij += fij[2*j];
        maij = mij - (double) fij[2*j]/2.0;
        bj = j+1;
        baj = bj - 0.5;
        
        if (j < l-1)
        {
            tmp = (double) nsum * mij - (double) Na * bj;
            innerSum = innerSum + tmp * tmp /
                (bj * ((double) nsum - bj));
        }

        tmp = (double) nsum * maij - (double) Na * baj;
        aInnerSum = aInnerSum + tmp * tmp /
            (baj * (nsum - baj) - nsum * 0.25);
    }

    AkN2 = AkN2 + innerSum / (Na * 1.0);
    AakN2 = AakN2 + aInnerSum / (Na * 1.0);

    //Start with array b
    mij = 0;
    maij = 0;
    innerSum = 0;
    aInnerSum = 0;

    for (int j=0; j<l; j++)
    {
        mij += fij[2*j+1];
        maij = mij - (double) fij[2*j+1]/2.0;
        bj = j+1;
        baj = bj - 0.5;

        if (j < l-1)
        {
            tmp = (double) nsum * mij - (double) Nb * bj;
            innerSum = innerSum + tmp * tmp /
                (bj * ((double) nsum - bj));
        }

        tmp = (double) nsum * maij - (double) Nb * baj;
        aInnerSum = aInnerSum + tmp * tmp /
            (baj * (nsum - baj) - nsum * 0.25);
    }

    AkN2 = AkN2 + innerSum / (Nb * 1.0);
    AakN2 = AakN2 + aInnerSum / (Nb * 1.0);


    AkN2 = AkN2 / (double)nsum;
    AakN2 = (nsum-1)*AakN2 / (nsum*nsum*1.0);

    double A2 = AkN2 - 1;
    A2 /= sigmaNorm;
   
//    std::cout << "Stat = " << A2  << "\n";
    double pvalue = PValueADKSamples(A2);

    return pvalue;
}
#endif //FRINGE_AD2UNIQUE_H
