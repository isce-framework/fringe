#include "filter.hpp"

int main(int argc, char **argv)
{
    int nwin=24;
    int npad = 8;
    float alpha = 0.2f;
    float lowpass = 30.0f;
    float beta = 0.5f;


    FilterWorker ant(nwin, npad, alpha, beta, lowpass);
    int nsize = ant.nfft();
    std::complex<float> cJ(0.0f, 1.0f);

    arma::cx_fmat indata( nwin, nwin);
    for(int ii=0; ii< nwin; ii ++)
    {
        for(int jj=0; jj < nwin; jj++)
        {

            indata(ii,jj) = std::exp(cJ * (0.0314f*ii + 0.02f*jj));
        }
    }

//    std::cout << "Original input matrix:  \n";
//    indata.submat(0,0,4,4).print();

//    std::cout << "Low pass filter: \n";
//    ant.LP.submat(0,0,4,4).print();


    ant.chip.submat(0, 0, nwin-1, nwin-1) = indata;
    ant.filter();

    arma::fmat angle = arma::arg(indata) - arma::arg(ant.chip.submat(0,0,nwin-1, nwin-1));

//    ant.chip.submat(0,0,4,4).print();

    angle.submat(0,0,4,4).print();


    std::cout << "Final coherence: " << arma::sum( arma::sum( arma::exp( cJ * arma::conv_to<arma::cx_fmat>::from(angle)))) / (1.0f * nwin * nwin) << "\n";

}
