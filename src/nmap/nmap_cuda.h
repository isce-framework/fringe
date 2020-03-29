//
// Author: Piyush Agram
// Copyright: 2018
//


#ifndef FRINGE_NMAP_CUDA_H
#define FRINGE_NMAP_CUDA_H

void lockGPU();
void unlockGPU();
void nmapProcessBlock(float *amp, unsigned char *msk,
                      int cols, int lines, int bands,
                      int *cnt, unsigned int *wmask,
                      int wtslen, double pval,
                      int Nx, int Ny);


#endif
