#include <iostream>
#include <complex>

/* This header file is just to include interface to LAPACK's
 * cheevr function to determine the largest eigen value of a 
 * Hermitian matrix - http://www.netlib.org/lapack/explore-html/df/d79/cheevr_8f_source.html */



/* Function signature
*  SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
*                          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
*                          RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBZ, RANGE, UPLO
*       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK,
*      $                   M, N
*       REAL               ABSTOL, VL, VU
*       ..
*       .. Array Arguments ..
*       INTEGER            ISUPPZ( * ), IWORK( * )
*       REAL               RWORK( * ), W( * )
*       COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * )
*/

//Confirmed signature against scipy/linalg/_lapack_subroutines.h
extern "C" 
{
    void zheevr_(char *JOBZ, char* RANGE, char* UPLO,
                    int *N, std::complex<double> *A, int *LDA,
                    double* VL, double* VU, int* IL, int* IU,
                    double *ABSTOL, int *M, double *W,
                    std::complex<double>* Z, int* LDZ, int* ISUPPZ,
                    std::complex<double>* WORK, int* LWORK,
                    double *RWORK, int *LRWORK,
                    int* IWORK, int* LIWORK, int* INFO);

    
/*    void zheev_(char *JOBZ, char* UPLO, int *N, std::complex<double> *A,
                    int *LDA, double *W, std::complex<double> *WORK,
                    int *LWORK, double *RWORK, int *INFO); */

    void zpotrf_(char *UPLO, int *N, std::complex<double> *A,
                    int *LDA, int *info);
    
    void zpotri_(char *UPLO, int *N, std::complex<double> *A, 
                    int* LDA, int *INFO);

/*    void dsyev_(char *JOBZ, char* UPLO, int* N,
                    double* A, int* LDA, double* W,
                    double* WORK, int* LWORK, int *INFO);

    void dsymm_(char *SIDE, char *UPLO, int *M, int *N,
                    double *A, int *LDA, double *B, int *LDB,
                    double *beta, double *C, double *LDC); */

}

struct EVWorker
{

    int size;           //Will store order of matrix
    double *eigval;     //To store eigen values
    std::complex<double> *eigvec;     //To store the eigen vectors

    int *isuppz;         //Worker array
    int Lisuppz;

    std::complex<double> *cxwork;    //Worker array
    int Lcxwork;

    double *rwork;        //Worker array
    int Lrwork;

    int  *iwork;      //Worker array
    int Liwork; 

    double *fwork;      //Worker array
    int Lfwork;

    EVWorker(){};
    EVWorker(int n);
    ~EVWorker();

    inline void prepare(int );
    inline int eigenWithIndex(std::complex<double>*, int, bool);
    inline int largestEigen(std::complex<double>*, bool);
    inline int smallestEigen(std::complex<double>*, bool);
    inline int positiveDefiniteInverse(std::complex<double>*);

};


EVWorker::EVWorker(int n)
{
    prepare(n);
}

EVWorker::~EVWorker()
{
    delete [] eigval;
    delete [] eigvec;
    delete [] isuppz;
    delete [] cxwork;
    delete [] rwork;
    delete [] iwork;
    delete [] fwork;
}

inline void EVWorker::prepare(int n)
{
    size = n;
    Lisuppz = 2*n;
    Lcxwork = 2*n;
    Lrwork = 24*n;
    Liwork = 10*n;
    Lfwork = n*n;

    eigval = new double[size];
    eigvec = new std::complex<double>[size*size];
    isuppz = new int[Lisuppz];
    cxwork = new std::complex<double>[Lcxwork];
    rwork = new double[Lrwork];
    iwork = new int[Liwork];
    fwork = new double[Lfwork];

}


inline int EVWorker::eigenWithIndex(std::complex<double> *A, int index, bool vecflag)
{
        
    char JOBZ;   //Compute Eigen Values and Vectors
    char RANGE = 'I';   //Compute limited set of solutions
    char UPLO = 'U';    //Compute assuming upper matrix contains data
    int N = size;       //Order of matrix

    int LDA = size;     //LDA is same as N
    double VL = 0.0;    //Not referenced
    double VU = 0.0;    //Not reference
    int IL = index;      //Get the specific eigen value.
    int IU = index;      //Get the specific eigen value.

    double abstol = 1.0e-6;  //Absolute tolerance

    int outM;       //Number of eigenvalues found

    int LDZ = size; //LDZ is same as N
    int info;       //For return status

    if (vecflag)    //If both values and vector are desired
    {
        JOBZ = 'V';
    }
    else            //If only values are needed
    {
        JOBZ = 'N';
    }

    /*std::cout << "Calling zheevr \n";
    std::cout << "JOBZ = " << JOBZ << "\n";
    std::cout << "RANGE = " << RANGE << "\n";
    std::cout << "UPLO = " << UPLO << "\n";
    std::cout << "N = " << N << "\n";
    std::cout << "A = " << A << "\n";
    std::cout << "VL = " << VL << "\n";
    std::cout << "VU = " << VU << "\n";
    std::cout << "IL = " << IL << "\n";
    std::cout << "IU = " << IU << "\n";
    std::cout << "ABSTOL = " << abstol << "\n";
    std::cout << "M = " << &outM << "\n";
    std::cout << "W = " << eigval << "\n";
    std::cout << "Z = " << eigvec << "\n";
    std::cout << "LDZ = " << LDZ << "\n";
    std::cout << "ISUPPZ = " << isuppz << "\n";
    std::cout << "(Lisuppz) = " << Lisuppz << "\n";
    std::cout << "WORK = " << cxwork << "\n";
    std::cout << "LWORK = " << Lcxwork << "\n";
    std::cout << "RWORK = " << rwork << "\n";
    std::cout << "LRWORK = " << Lrwork << "\n";
    std::cout << "IWORK = " << iwork << "\n";
    std::cout << "LIWORK = " << Liwork << "\n";*/

    zheevr_(&JOBZ, &RANGE, &UPLO, &N, A,
            &LDA, &VL, &VU, &IL, &IU,
            &abstol, &outM, eigval, eigvec,
            &LDZ, isuppz, cxwork, &Lcxwork,
            rwork, &Lrwork, iwork, &Liwork,
            &info); 

    /*std::cout << "INFO = " << info << "\n";
    std::cout << "M = " << outM << "\n";
    std::cout << "val = " << eigval[0] << "\n";
    std::cout << eigvec[0] << " " << eigvec[1] << " " << eigvec[2] << "\n";

    std::cout << "cxwork = " << Lcxwork << "  " << cxwork[0] << "\n"
              << "rwork  = " << Lrwork << "  " << rwork[0] << "\n"
              << "iwork  = " << Liwork << "  " << iwork[0] << "\n";*/

    return info;
}


inline int EVWorker::largestEigen(std::complex<double> *A, bool flag)
{
    int index=size;
    return eigenWithIndex(A, index, flag);
}

inline int EVWorker::smallestEigen(std::complex<double> *A, bool flag)
{
    int index = 1;
    return eigenWithIndex(A, index, flag);
}

inline int EVWorker::positiveDefiniteInverse(std::complex<double> *A)
{
    int info;

    {
        char UPLO='U';
        int N = size;
        int LDA=size;

        zpotrf_(&UPLO, &N, A,
                &LDA, &info);
    }

    if (info !=0) return info;

    {
        char UPLO = 'U';
        int N = size;
        int LDA = size;
        zpotri_(&UPLO, &N, A,
               &LDA, &info);
    }

    //Make the matrix hermitian
    for(int ii=0; ii<size; ii++)
    {
        for(int jj=ii+1; jj<size; jj++)
        {
            std::complex<double> val = A[jj*size + ii];
            A[ii*size + jj] = std::conj(val);
        }
    }

    return info;
}
