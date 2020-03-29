#ifndef FRINGE_CUDA_ULONGMASK_H
#define FRINGE_CUDA_ULONGMASK_H

#include <cuda_runtime.h>


__device__ inline void setBit(unsigned int *arr, 
                        int ii, int jj, 
                        int Nx, int Ny, 
                        bool flag)
{
    //Unity
    unsigned int U1 = 1;

    //Flat index into 2D array
    int bit = (ii+Ny)*(2*Nx+1)+jj+Nx;

    //Byte number into arr to modify
    int num = bit/32;

    //Bit number of specific byte to modify
    bit = bit%32;

    if (flag)
    {
        arr[num] |= ((U1)<<bit);
    }
    else
    {
        arr[num] &= (~((U1)<<bit));
    }
}

__device__ inline bool getBit(unsigned int *arr,
                        int ii, int jj,
                        int Nx, int Ny)
{
    //Unity
    unsigned int U1 = 1;

    //Flat index into 2D array
    int bit = (ii+Ny)*(2*Nx+1)+jj+Nx;

    //Byte number into arr to modify
    int num = bit/32;

    //Bit number of specific byte to modify
    bit = bit%32;

    return (U1 == ((arr[num] >> bit) & U1));

}
#endif
