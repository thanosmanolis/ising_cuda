#ifndef CUDA_H
#define CUDA_H

#include "ising.h"

/*
*************************************************
*    @file   cuda.h                             *
*    @author amanolis <amanolis@ece.auth.gr>    *
*    @date   Sat Jan 10 16:03:30 2020           *
*    @brief  Cuda functions                     *
*************************************************
*/

/*
************************************************************************************
*    CUDA Kernel                                                                   *
*                                                                                  *
*    - param n            Number of lattice points per dim            [scalar]     *
*    - param gpu_w        GPU array that stores the weight matrix     [5-by-5]     *
*    - param gpu_G        GPU array that stores the atomic spins      [n-by-n]     *
*    - param gpu_G_new    GPU array that stores the updated values    [n-by-n]     *
*                                                                                  *
*    NOTE: All matrices are stored in row-major format                             *
************************************************************************************
*/

__global__ void kernel(int n, double* gpu_w, int* gpu_G, int* gpu_G_new);

#endif
