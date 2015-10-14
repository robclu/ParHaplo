// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_GPU_H
#define PARHAPLO_TREE_GPU_H

#include "read_info_gpu.cuh"
#include "snp_info_gpu.h"

namespace haplo {
    
class TreeGpu {
public:
    //-------------------------------------------------------------------------------------------------------
    using binary_vector         = BinaryVector<2>;
    using read_info_type        = ReadInfo;
    using read_info_container   = thrust::host_vector<read_info_type>;
    using small_type            = uint8_t;
    //-------------------------------------------------------------------------------------------------------
private:
    binary_vector&              _data;
    read_info_container&        _read_info;
    
    // ------------------------------------------ DEVICE ----------------------------------------------------
    
    small_type*                 _device_data;
    read_info_type*             _device_read_info;
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    TreeGpu(binary_vector& data, read_info_container& read_info);
    
    CUDA_H
    ~TreeGpu() 
    {
        cudaFree(_device_data);
    }
};


// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

TreeGpu::TreeGpu(binary_vector& data, read_info_container& read_info)
: _data(data), _read_info(read_info)  
{
    // Allocate memory on the device for the data
    cudaMalloc((void**)&_device_data, sizeof(small_type) * _data.size());
//    cudaMemcpy(_device_data, _data, sizeof(binary_vector_big) * _data.size(), cudaMemcpyHostToDevice); 
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


