// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo cuda definitions
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_CUDA_DEFS_H
#define PARHAPLO_CUDA_DEFS_H

#ifdef __CUDACC__
    #define CUDA_HD __host__ __device__
    #define CUDA_H  __host__
    #define CUDA_D  __device__
    #define SHARED  __shared__
    #define ALIGN(x) __align__(x)
#else
    #define CUDA_HD
    #define CUDA_H
    #define CUDA_D
    #define SHARED
    #define ALIGN(x) 
#endif

#endif      // PARAHAPLO_CUDA_DEFS_H
