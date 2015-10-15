// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_GPU_H
#define PARHAPLO_TREE_GPU_H

#include "tree_kernels_gpu.cu"
#include "tree_node_manager.h"

namespace haplo {
    
class TreeGpu {
public:
    //-------------------------------------------------------------------------------------------------------
    using binary_vector         = BinaryVector<2>;
    using read_info_type        = ReadInfo;
    using read_info_container   = thrust::host_vector<read_info_type>;
    using snp_info_type         = SnpInfoGpu;
    using snp_info_container    = thrust::host_vector<snp_info_type>;
    using small_type            = uint8_t;
    using node_container        = NodeManagerGpu;
    //-------------------------------------------------------------------------------------------------------
private:
    binary_vector&              _data;
    read_info_container&        _read_info;
    snp_info_container          _snp_info;
    
    // ------------------------------------------ DEVICE ----------------------------------------------------
    
    small_type*                 _device_data;
    read_info_type*             _device_read_info;
    snp_info_type*              _device_snp_info;
    node_container              _device_nodes;
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    TreeGpu(binary_vector& data, read_info_container& read_info, snp_info_container, const size_t nodes);
    
    CUDA_H
    ~TreeGpu() 
    {
        cudaFree(_device_data);
        cudaFree(_device_read_info);
        cudaFree(_device_snp_info);
        cudaFree(_device_nodes._nodes);
    }
};


// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

TreeGpu::TreeGpu(binary_vector& data, read_info_container& read_info, snp_info_container snp_info,
                 const size_t nodes)
: _data(data), _read_info(read_info), _snp_info(snp_info), _device_nodes(nodes)
{
    cudaError_t error;
    
    // Copy data to the device
    thrust::host_vector<small_type> data_gpu_format = data.to_binary_vector();
    cudaMalloc((void**)&_device_data, sizeof(small_type) * data.size());
    cudaMemcpy(_device_data, thrust::raw_pointer_cast(&data_gpu_format[0]), 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice); 

    // Copy the read data to the device 
    cudaMalloc((void**)&_device_read_info, sizeof(read_info_type) * read_info.size());
    cudaMemcpy(_device_read_info, thrust::raw_pointer_cast(&_read_info[0]),
                    sizeof(read_info_type) * read_info.size(), cudaMemcpyHostToDevice); 
    
    // Copy SNP data to the device 
    cudaMalloc((void**)&_device_snp_info, sizeof(snp_info_type) * snp_info.size());
    cudaMemcpy(_device_snp_info, thrust::raw_pointer_cast(&_snp_info[0]), 
                    sizeof(snp_info_type) * snp_info.size(), cudaMemcpyHostToDevice); 
   
    size_t* result; 
    cudaMalloc((void**)&result, sizeof(size_t));
    
    // Try and invoke the kernel
    search_tree<<<128,1>>>(_device_data, data.size(), _device_read_info, read_info.size(),
                           _device_snp_info, snp_info.size(), result);
   
    size_t host_result;
    cudaMemcpy(&host_result, result, sizeof(size_t), cudaMemcpyDeviceToHost); 
    
    std::cout << host_result << "\n\n";
    
    // Allocate space on the device for the tree nodes
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    
    std::cout << "Free Memory: " << free_mem << "\n";
    size_t num_nodes = 0.7* free_mem / sizeof(TreeNode);
    _device_nodes.num_nodes() = num_nodes;
    error = cudaMalloc((void**)&_device_nodes._nodes, sizeof(TreeNode) * num_nodes);
    
    if (error != cudaSuccess) std::cout << "Error allocating memory\n";
    
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "Free Memory: " << free_mem << "\n";

 
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


