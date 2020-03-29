//
// Author: Piyush Agram
// Copyright 2018

#include "fringe/cuda/cudaUtils.h"
#include "fringe/cuda/ulongmask.h"
#include "KS2sample_cuda.h"
#include <math.h>
#include <iostream>

#define THRD_PER_BLOCK 96

//Constant memory for constant input values
//Done this way in topozero. Need to understand
//why this cannot be put in part of the struct
__constant__ double d_inpts_pval[1];
__constant__ int d_inpts_int[6];
//0: numCols
//1: numLines
//2: nBands
//3: Nx
//4: Ny
//5: wtslen

/*******Advance function declarations***********/
__global__ void runSortAmp(float *amp);
__global__ void findNeighbors(const float *amp, const unsigned char *mask,
                            unsigned int *wts, int * count);

/**********Inplace sorting algorithm************/
__device__ inline void bruteForceSort(float *arr)
{
    int ii,jj;
    float temp;

    for(jj=0; jj < (d_inpts_int[2]-1); jj++)
    {
        for (ii=0; ii < (d_inpts_int[2]-jj-1); ii++)
        {
            temp = fminf(arr[ii], arr[ii+1]);
            arr[ii+1] = fmaxf(arr[ii], arr[ii+1]);
            arr[ii] = temp;
        }
    }
};

/***** Structure for all GPU/Host handling data *****/
struct gpuParams
{
    //These are meant to be inputs
    //Pointers refer to pointers on device
    float *amplitude; //Flattened array of amplitude values
    unsigned char *mask; //Flattened mask array
    
    //These are host end variables
    int numCols;      //Number of cols in a line
    int numLines;     //Number of lines in a block
    int nBands;       //Number of bands
    int numWtsBands;    //Number of uint32 bands for nmap
    int Nx;           //Half window width in X 
    int Ny;           //Half window width in Y

    //These are meant to be outputs
    //Pointers refer to pointers on device
    int *count;         //Count of number of neighbors
    unsigned int *wts;  //Weights

    //Constructor and destructor
    gpuParams(int cols, int lines, int bands,
              int nx, int ny, int nwts);
    ~gpuParams();

    //Methods to help
    void allocateArrays();
    void setConstants();
    void deallocateArrays();
    void setInputs(float *, unsigned char *);
    void getOutputs(int *, unsigned int *);
    void sortAmplitude();
    void process(double pval);
};

gpuParams::gpuParams(int cols, int lines, int bands,
                     int nx, int ny, int nwts):
                     numCols(cols),numLines(lines),
                     nBands(bands),
                     Nx(nx), Ny(ny),
                     numWtsBands(nwts)
{
    //Ensure the memory on GPU is allocated for this
    allocateArrays();

    //Assign constants to the global arrays
    setConstants();
}

gpuParams::~gpuParams()
{
    //Ensure memory on GPU is released
    deallocateArrays();
}

//Allocate memory on the GPU
void gpuParams::allocateArrays()
{

    size_t nPix = numCols * numLines;

    //Allocate memory for input amplitude
    gpuErrChk ( cudaMalloc((float**)&amplitude, (nBands*nPix)*sizeof(float)));
   gpuErrChk( cudaMemset(amplitude, 0, (nBands*nPix)*sizeof(float)));

    //Allocate memory for input mask
    gpuErrChk( cudaMalloc((unsigned char**)&mask, (nPix)*sizeof(unsigned char)));
    gpuErrChk( cudaMemset(mask, 0, nPix*sizeof(unsigned char)));

    //Allocate memory for output count
    gpuErrChk( cudaMalloc((int**)&count, (nPix)*sizeof(int)));
    gpuErrChk( cudaMemset(count, 0, nPix*sizeof(int)));

    //Allocate memory for weights
    gpuErrChk( cudaMalloc((unsigned int**)&wts, (nPix*numWtsBands)*sizeof(unsigned int)));
   gpuErrChk( cudaMemset(wts, 0, nPix*numWtsBands*sizeof(unsigned int))  );

}

void gpuParams::setConstants()
{
    int constants[6];
    constants[0] = numCols;
    constants[1] = numLines;
    constants[2] = nBands;
    constants[3] = Nx;
    constants[4] = Ny;
    constants[5] = numWtsBands;

    gpuErrChk( cudaMemcpyToSymbol(d_inpts_int, constants, (6*sizeof(int))));

    //int readback[6];
    //gpuErrChk( cudaMemcpyFromSymbol(readback, d_inpts_int, (6*sizeof(int))));
    //std::cout << "Ncols = " << readback[0] << "\n"
    //          << "Nlines = " << readback[1] << "\n"
    //          << "Nbands = " << readback[2] << "\n"
    //          << "Nx = " << readback[3] << "\n"
    //          << "Ny = " << readback[4] << "\n"
    //          << "Nwts = " << readback[5] << "\n";

}

//Deallocate memory on GPU
void gpuParams::deallocateArrays()
{
    //Free amplitude
    gpuErrChk( cudaFree(amplitude));

    //Free mask
    gpuErrChk( cudaFree(mask));

    //Free output
    gpuErrChk( cudaFree(count));

    //Free weights
    gpuErrChk( cudaFree(wts));
}

//Pass amplitude and mask to GPU
void gpuParams::setInputs(float *amp,
                          unsigned char *msk)
{
    size_t nPix = numCols * numLines;
    
    //Copy amplitude to GPU
    gpuErrChk( cudaMemcpy(amplitude, amp, (nPix*nBands*sizeof(float)),
            cudaMemcpyHostToDevice));

    //Copy mask to GPU
    gpuErrChk( cudaMemcpy(mask, msk, (nPix*sizeof(unsigned char)),
            cudaMemcpyHostToDevice));
}

//Get count and wts from GPU
void gpuParams::getOutputs(int *cnt,
                           unsigned int *wmask)
{
    size_t nPix = numCols * numLines;

    //Copy count to host
    gpuErrChk( cudaMemcpy(cnt, count, (nPix*sizeof(int)),
            cudaMemcpyDeviceToHost));

    //copy wts to host
    gpuErrChk( cudaMemcpy(wmask, wts, (nPix*numWtsBands*sizeof(unsigned int)), cudaMemcpyDeviceToHost));

}

//Sort amplitudes
void gpuParams::sortAmplitude()
{
    int numPix = numCols * numLines;
    dim3 block(THRD_PER_BLOCK);
    dim3 grid((numPix + (THRD_PER_BLOCK-1))/THRD_PER_BLOCK);

    /*if ((grid.x * THRD_PER_BLOCK) > numPix)
    {
        std::cout << " Number of empty threads = " << ((grid.x * THRD_PER_BLOCK) - numPix) << "\n";
    }*/

    runSortAmp <<<grid, block>>>(amplitude);

    //Track errors and synchronize
    gpuErrChk( cudaGetLastError());
    gpuErrChk( cudaDeviceSynchronize());
}

//Find neighbors
void gpuParams::process(double pval)
{
    //Copy the threshold value to device
    gpuErrChk( cudaMemcpyToSymbol(d_inpts_pval, &pval, sizeof(double)));

    int numPix = numCols * numLines;
    dim3 block(THRD_PER_BLOCK);
    dim3 grid((numPix + (THRD_PER_BLOCK-1))/THRD_PER_BLOCK);

    /*if ((grid.x * THRD_PER_BLOCK) > numPix)
    {
        std::cout << "Number of empty threads = " << ((grid.x * THRD_PER_BLOCK) - numPix) << "\n";
    }*/

    findNeighbors <<<grid, block>>>(amplitude, mask,
                                    wts, count);

    //Track errors and synchronize
    gpuErrChk(cudaGetLastError());
    gpuErrChk(cudaDeviceSynchronize());

}

/*****End of structure************/

/********** Actual Kernel function *************/
//This method is to sort a single pixel 
__global__ void runSortAmp(float *amp)
{
    //Pixel number 
    int pixel = (blockDim.x * blockIdx.x) + threadIdx.x;

    //Make sure count is within limits
    //i.e, pixel < numCols * numLines
    if (pixel < (d_inpts_int[0] * d_inpts_int[1]))
    {
        //Offset to pixel = pixel * nbands
        bruteForceSort(amp+(pixel*d_inpts_int[2])); 
    }
}


//This method is to identify neighbors for a single pixel
__global__ void findNeighbors(const float *amp, const unsigned char *mask,
                    unsigned int *wts, int *count)
{
    //Temporary variables needed
    int refii, refjj;
    int qq,ii,jj;
    double prob;
    const float *refpix;
    const float *cenpix;
    unsigned int *weight;

    //Pixel number
    int pp = (blockDim.x * blockIdx.x) + threadIdx.x;

    //Make sure count is within limits
    //i.e, pixel < numCols * numLines
    if (pp < (d_inpts_int[0] * d_inpts_int[1]))
    {
        if( mask[pp] != 0)
        {
            cenpix = amp + (d_inpts_int[2] * pp);
            weight = wts + (d_inpts_int[5] * pp);
            
            for(ii=-d_inpts_int[4]; ii<=d_inpts_int[4]; ii++)
            {
                refii = (pp/d_inpts_int[0]) + ii;

                for (jj=-d_inpts_int[3]; jj<=d_inpts_int[3]; jj++)
                {
                    refjj = (pp%d_inpts_int[0]) + jj;

                    if ((refii < d_inpts_int[1]) && (refii >=0) && 
                        (refjj < d_inpts_int[0]) && (refjj >=0)) 
                    {
                        qq = refii * d_inpts_int[0] + refjj;
                        refpix = amp + (d_inpts_int[2] * qq); 
                            
                        if (mask[qq] != 0)
                        {
                            //Count same pix as neighbor
                            if (pp == qq)
                            {
                                count[pp] += 1;
                                setBit(weight,
                                        0, 0,
                                        d_inpts_int[3],
                                        d_inpts_int[4],
                                        true);
                            }
                            else
                            {
                                prob = KS2test(cenpix, refpix,
                                                d_inpts_int[2]);
                                   
                                if (prob >= d_inpts_pval[0])
                                {
                                    count[pp] += 1;

                                    setBit( weight,
                                            ii,jj,
                                            d_inpts_int[3],
                                            d_inpts_int[4],
                                            true);
                                } //if prob > thresh
                            }  //if not same pixel
                        } //if ref pixel is not masked
                    } //if ref pixel is within limits
                } //loop over jj
            } //loop over ii 
            /*count[pp] = 0;
            for(ii=0; ii< d_inpts_int[2]; ii++)
                count[pp] += (cenpix[ii] == 0);*/

        } //if pixel is not masked
    } //if pixel is within limits
}

                
/********* End of actual kernel function **********/


//Actual interface to nmap.cpp
//This is the only function that is used by parent code directly.
void nmapProcessBlock(float *amp, unsigned char *msk,
                      int cols, int lines, int bands,
                      int *cnt, unsigned int *wmask,
                      int wtslen, double pval,
                      int Nx, int Ny)
{

    //Create structure to handle interaction with GPU
    struct gpuParams pars( cols,lines, bands,
                           Nx, Ny, wtslen);

    //Copy inputs to GPU
    pars.setInputs(amp, msk);

    //Sort amplitudes since stats is between histograms
    pars.sortAmplitude();

    //Process the block
    pars.process(pval);

    //Get outputs from GPU
    pars.getOutputs(cnt, wmask);

}

//Wrappers for GPU access
void lockGPU()
{
    getGPUDevice(0);
}

void unlockGPU()
{
    releaseGPUDevice();
}
