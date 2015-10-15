// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class kernels
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_KERNELS_INTERFACE_GPU_H
#define PARHAPLO_TREE_KERNELS_INTERFACE_GPU_H

#include "tree_kernels_gpu.cu"

namespace haplo {

template <size_t ThreadsX = 128, size_t ThreadsY = 1, size_t ThreadsZ = 1> 
size_t search_tree_interface(uint8_t* data       , size_t data_size, 
                             ReadInfo* read_info , size_t read_size, 
                             SnpInfoGpu* snp_info, size_t snps_size)
{
    size_t* result; 
    cudaMalloc((void**)&result, sizeof(size_t));
    
    // Invoke the kernel
    search_tree<<<ThreadsX,ThreadsY>>>(_device_data, data.size(), _device_read_info, read_info.size(),
                                     _device_snp_info, snp_info.size(), result);
   
    size_t host_result;
    cudaMemcpy(&host_result, result, sizeof(size_t), cudaMemcpyDeviceToHost); 
    
    return host_result;
}


};              // End namespace haplo
#endif          // PARAHAPLO_TREE_KERNELS_GPU_H
    