#include <armadillo>
#include <complex>
#include <fftw3.h>

struct FilterWorker
{
    //Members
    const int nwin;       //Window size for data
    const int npad;       //Zero padding for the FFTs

    //Filter parameters
    const float alpha;      //Power of spectral filter
    const float beta;       //Strength of spectral filter
    const float lowpass;    //Wavelength of low pass filter


    //Private arrays
    arma::cx_fmat chip;      //Chip of data
    arma::cx_fmat SP;        //Spectral filter matrix
    arma::cx_fmat LP;        //Low pass filter matrix
    arma::cx_fmat  GW;        //Gaussian smoothing window
    arma::cx_fmat work;       //Temporary arrays
    arma::cx_fmat workpad;    


    //Private plans
    fftwf_plan fwdplan;         //Stores forward plan
    fftwf_plan invplan;         //Stores inverse plan

    fftwf_plan fwdpadplan;      //FFT with padding
    fftwf_plan invpadplan;      //IFFT with padding

    //Functions
    inline int nfft()     //Sum of data window + padding
    {
        return npad + nwin;
    };

    inline int nfftpad()    //Sum of data window + padding + gaussian padding
    {
        return nfft() + 8;
    };

    inline void setupPlans();     //Set up FFT plans
    inline void destroyPlans();   //Destroy FFT plans
    inline void setupLowPass();   //Setup Butterwortjh filter
    inline void setupGaussian();  //Setup Gaussian window
    inline void filter();         //Do the actual filtering

    inline void zeros()
    {
        SP.zeros();
        chip.zeros();
    };

    //Constructor
    FilterWorker(int nwin_,int npad_, float a_, float b_, float l_):
        nwin(nwin_),
        npad(npad_),
        alpha(a_),
        beta(b_),
        lowpass(l_),
        chip(npad_ + nwin_, npad_ + nwin_),
        SP(npad_ + nwin_, npad_ + nwin_),
        LP(npad_ + nwin_, npad_ + nwin_),
        GW(npad_ + nwin_ + 8, npad_ + nwin_ + 8),
        work(npad_ + nwin_, npad_ + nwin_),
        workpad(npad_ + nwin_ + 8, npad_ + nwin_ + 8)
    {
        setupPlans();
        setupLowPass();
        setupGaussian();
    };

    //Destructor
    ~FilterWorker()
    {
        destroyPlans();
    };
};

//In memory setting of imaginary part to zero
inline void setImagPartToZero(arma::cx_fmat &mat)
{
    arma::fmat spread  = arma::fmat( (float*) mat.memptr(), 2 * mat.n_rows,
                                        mat.n_cols, false);
    spread.rows(arma::regspace<arma::uvec>(1,2,spread.n_rows-1)).zeros();
}


//See if any element is non-zero
inline bool hasNonZero(arma::cx_fmat &mat)
{
    arma::fvec spread = arma::fvec( (float*) mat.memptr(), 2*mat.n_rows*mat.n_cols, false);
    return arma::any(spread);
}


//Find median of complex matrix
inline std::complex<float> findMedian(arma::cx_fmat &mat)
{
    arma::cx_fvec spread = arma::cx_fvec( (std::complex<float>*) mat.memptr(), mat.n_rows*mat.n_cols, false);
    return arma::median(spread);
}

//FFTshift to handline padding. Not to confuse with actual FFTshift.
inline void FFTrepack(arma::cx_fmat &src, arma::cx_fmat &dst, int nsize)
{
    dst.zeros();
    int halfsize = nsize/2;

    //Quad BR to TL
    dst.submat(0,0, halfsize - 1, halfsize-1) = src.submat(halfsize, halfsize, nsize-1, nsize-1);
    
    //Quad TL to BR
    dst.submat(halfsize, halfsize, nsize-1, nsize-1) = src.submat(0,0,halfsize-1, halfsize-1);

    //Quad TR to BL
    dst.submat(halfsize, 0, nsize-1, halfsize-1) = src.submat(0, halfsize, halfsize-1, nsize-1);

    //Quad BL to TR
    dst.submat(0, halfsize, halfsize-1, nsize-1) = src.submat(halfsize, 0, nsize-1, halfsize-1);
}

//Create in place FFT plans
inline void FilterWorker::setupPlans()
{
    int nsize = nfft();

    int nlarge = nfftpad();

    //Forward FFT plan
    fwdplan = fftwf_plan_dft_2d(nsize, nsize, (fftwf_complex*) chip.memptr(), (fftwf_complex*) chip.memptr(), FFTW_FORWARD, FFTW_MEASURE);

    //Inverse FFT plan
    invplan = fftwf_plan_dft_2d(nsize, nsize, (fftwf_complex*) chip.memptr(), (fftwf_complex*) chip.memptr(), FFTW_BACKWARD, FFTW_MEASURE);

    //Forward FFT plan for gaussian filtering
    fwdpadplan = fftwf_plan_dft_2d(nlarge, nlarge, (fftwf_complex*) workpad.memptr(), (fftwf_complex*) workpad.memptr(), FFTW_FORWARD, FFTW_MEASURE);

    //Inverse FFT plan
    invpadplan = fftwf_plan_dft_2d(nlarge, nlarge, (fftwf_complex*) workpad.memptr(), (fftwf_complex*) workpad.memptr(), FFTW_BACKWARD, FFTW_MEASURE);


}


//Destroy FFT plans
inline void FilterWorker::destroyPlans()
{
    fftwf_destroy_plan(fwdplan);
    fftwf_destroy_plan(invplan);
    fftwf_destroy_plan(fwdpadplan);
    fftwf_destroy_plan(invpadplan);
}


inline void FilterWorker::setupLowPass()
{
    int nsize = nfft();
    
    if (lowpass <= 0.0f)
    {
        LP.zeros();
    }
    else
    {
        LP.ones();
        float freq0 = 1.0f / lowpass;
        for(int ii=0; ii< nsize; ii++)
        {
            float ff = (ii >= (nsize/2))? (ii-nsize) : ii;
            ff /= (1.0f*nsize);

            float butter = 1.0f / (1.0f + std::pow(ff/freq0, 10.0f));
            
//            std::cout << "Butter: " << ff << " " << std::pow(ff/lowpass, 10.0f) << "\n";
            LP.row(ii) *= std::complex<float>(butter,0.0f);
            LP.col(ii) *= std::complex<float>(butter,0.0f);
        }
    }
}

inline void FilterWorker::setupGaussian()
{
    //gausswin(7)
    arma::fvec gwin = {0.04393693f,0.24935221f,0.70664828f,1.0f,0.70664828f,0.24935221f,0.04393693f};

    int nsize = nfftpad();
    arma::cx_fvec padded(nsize);
    padded.zeros();

    //Wrap around
    for (int ii=0; ii<4; ii++)
        padded[ii].real(gwin[3+ii]);

    for (int ii=0; ii<3;ii++)
        padded[nsize-3+ii].real(gwin[ii]);


    //Symmetric 2D from 1D array
    workpad.ones();
    for (int ii=0; ii<nsize; ii++)
    {
        workpad.row(ii) *= padded[ii];
        workpad.col(ii) *= padded[ii];
    }

    //chip.submat(0,0,3,3).print();

    //Convert to spectra and store
    fftwf_execute(fwdpadplan);
    GW = workpad/(1.0f * nsize * nsize);
    workpad.zeros();

}


//Assumes that data has already been copied into the chip
inline void FilterWorker::filter()
{
    int nsize = nfft();
    int nlarge = nfftpad();


    //Zero padding
    chip.tail_cols(npad).zeros();
    chip.tail_rows(npad).zeros();

    //If chip is empty, skip all processing
    bool status = hasNonZero(chip);
    if (!status)
    {
        chip.zeros();
        return;
    }

    //If spectral component is desired.
    if (beta > 0.0f)
    {

        //Execute FFT
        fftwf_execute(fwdplan);

        //Copy to temporary array
        work.zeros();
        work.set_real(arma::abs(chip));

        //Execute FFTshift
        FFTrepack(work, workpad, work.n_rows);

        //convert abs value of spectra to time domain
        fftwf_execute(fwdpadplan);

        //multiply with gaussian spectra 
        workpad %= GW;

        //Inverse FFT to complete convolution
        fftwf_execute(invpadplan);

        //Execute IFFTshift
        FFTrepack(workpad, work, work.n_rows);


        //Setting imaginary part to zero
        setImagPartToZero(work);

        //Compute median
        std::complex<float> median = findMedian(work);
//        std::cout << "Median = " << median << "\n";

        //Scale by median if needed
        if (median != 0.0f)
        {
            work /= median;
        }

        //Raise to power
        work = arma::pow(work, alpha);
        work -= std::complex<float>(1.0f,0.0f);


        for(int ii=0; ii< work.n_elem; ii++)
        {
            float val = work[ii].real();
            val = (val > 0.0f)? val: 0.0f;
            work[ii].real(val);
        }

        work *= beta;
    }
    else
    {
        work.zeros();
    }

    //Combine low pass and spectral component here
    chip %= (work + LP)/(1.0f * nsize * nsize);

    /*std::cout << "Start \n";
    chip.submat(0,0,4,4).print();
    std::cout << "Mid \n";
    chip.submat(13,13,17,17).print();
    std::cout << "End \n";
    chip.submat(27,27,31,31).print();*/


    //Perform inverse FFT
    fftwf_execute(invplan);

    return;
}
 
