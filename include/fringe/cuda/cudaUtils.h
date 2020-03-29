//
//Author: Piyush Agram
//Copyright 2016
//

#ifndef FRINGE_CUDAUTILS_H
#define FRINGE_CUDAUTILS_H

#include <cuda_runtime.h>
#include <iostream>

//CUDA error handler
#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
       std::cerr<< "GPUassert: " << cudaGetErrorString(code) << " " << file << "  " << line << "\n";
      if (abort) exit(code);
   }
}

//Number of CUDA devices
int getNumDevices()
{
    int nDevices;
    gpuErrChk( cudaGetDeviceCount(&nDevices));
    
    return nDevices;
}


//Get GPU device
void getGPUDevice(int num=0)
{
    int ngpus = getNumDevices();

    if (num < ngpus)
    {
        //Get access to specific device
        gpuErrChk( cudaSetDevice(num));
    }
    else
    {
        std::cout << "Device index " << num << "outside range [0.."
                  << ngpus-1 << "] \n";
        exit(1);
    }
}

//Release GPU device
void releaseGPUDevice()
{
    gpuErrChk( cudaDeviceReset());
}


//Function to determine memory available on device
size_t getDeviceMem()
{
    size_t freeByte, totalByte;
    gpuErrChk( cudaMemGetInfo(&freeByte, &totalByte));
    totalByte = (totalByte/1e9) * 1e9; //Round down to nearest GB
    return totalByte;
}


#endif
