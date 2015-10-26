// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRAPH_GPU_H
#define PARHAPLO_GRAPH_GPU_H

#include "cuda_error.h"
#include "data.h"
#include "devices.hpp"
#include "graph.h"
#include "graph_kernels_gpu.cu"
#include "thrust/sequence.h"

#define EDGE_MEM_PERCENT 0.6f

namespace haplo {

// Specialization for gpu
template <typename SubBlockType>    
class Graph<SubBlockType, devices::gpu> {
public:
    //-------------------------------------------------------------------------------------------------------
    using binary_vector                 = BinaryVector<2>;
    using data_type                     = Data;
    using graph_type                    = internal::Graph;
    using read_info_type                = typename data_type::read_info_type;
    using read_info_container           = thrust::host_vector<read_info_type>;
    using snp_info_type                 = typename data_type::snp_info_type;
    using snp_info_container            = thrust::host_vector<snp_info_type>;
    using small_type                    = typename data_type::small_type;
    using small_container               = thrust::host_vector<small_type>;
    //-------------------------------------------------------------------------------------------------------
private:
    SubBlockType&               _sub_block;
    // -------------------------------------------- HOST ----------------------------------------------------
    binary_vector&              _data_cpu;
    read_info_container&        _read_info;
    snp_info_container          _snp_info;
    small_container             _haplotype;
    small_container             _alignments;
    size_t                      _snps;
    size_t                      _reads;
    size_t                      _device;
    size_t                      _nih_cols;
    
    // ------------------------------------------ DEVICE ----------------------------------------------------
    data_type                   _data_gpu; 
    graph_type                  _graph;
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    Graph(SubBlockType& sub_block, const size_t device);

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Destuctor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    ~Graph() 
    {
        cudaFree(_data_gpu.data     );
        cudaFree(_data_gpu.read_info);
        cudaFree(_data_gpu.snp_info );
        cudaFree(_graph.edges       );
    }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Solves the graph for the haplotypes
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void search();
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Sorts the edges
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void sort_edges(dim3 grid_size, dim3 block_size);
};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

template <typename SubBlockType>
Graph<SubBlockType, devices::gpu>::Graph(SubBlockType& sub_block, const size_t device)
: _sub_block(sub_block)                     , _data_cpu(sub_block.data())               , 
  _read_info(sub_block.read_info())         , _snp_info(sub_block.snp_info())           , 
  _haplotype(sub_block.snp_info().size())   , _alignments(sub_block.read_info().size()) , 
  _snps(sub_block.snp_info().size())        , _reads(sub_block.read_info().size())      ,
  _device(device)                           , _nih_cols(sub_block.nih_columns())        ,
  _data_gpu(sub_block.snp_info().size()     , sub_block.read_info().size())             
{
    cudaError_t error;
   
    std::cout << _data_gpu.reads << " " << _data_gpu.snps << "\n";
    
    // Copy the actual data to the device
    thrust::host_vector<small_type> data_gpu_format = _data_cpu.to_binary_vector();
    CudaSafeCall( cudaMalloc((void**)&_data_gpu.data, sizeof(small_type) * _data_cpu.size()) );
    CudaSafeCall( cudaMemcpy(_data_gpu.data, thrust::raw_pointer_cast(&data_gpu_format[0])   , 
                    sizeof(small_type) * _data_cpu.size(), cudaMemcpyHostToDevice)           );

    // Copy the read data to the device 
    CudaSafeCall( cudaMalloc((void**)&_data_gpu.read_info, sizeof(read_info_type) * _reads) );
    CudaSafeCall( cudaMemcpy(_data_gpu.read_info, thrust::raw_pointer_cast(&_read_info[0])  ,
                    sizeof(read_info_type) * _reads, cudaMemcpyHostToDevice)                );
    
    // Copy SNP data to the device 
    CudaSafeCall( cudaMalloc((void**)&_data_gpu.snp_info, sizeof(snp_info_type) * _snps)    );
    CudaSafeCall( cudaMemcpy(_data_gpu.snp_info, thrust::raw_pointer_cast(&_snp_info[0])    , 
                    sizeof(snp_info_type) * _snps, cudaMemcpyHostToDevice)                  );
 
    // Check that there is enough space for the edges
    size_t free_memory, total_memory;
    cudaMemGetInfo(&free_memory, &total_memory);
    
    // The number of total edges required
    const size_t num_edges = _reads * (_reads - 1) / 2;
    
    if (num_edges * sizeof(Edge) < free_memory * EDGE_MEM_PERCENT) { 
        // Allocate space for the edges of the graph
        CudaSafeCall( cudaMalloc((void**)&_graph.edges, sizeof(Edge) * num_edges) );
    } else { // Input is too big
        std::cerr << "Error : Input is too big for GPU =(\n";
    }
}

template <typename SubBlockType>
void Graph<SubBlockType, devices::gpu>::search()
{
    // The number of total edges required
    const size_t num_edges = _reads * (_reads - 1) / 2; 
    
    // Properties of the device
    int device; cudaDeviceProp device_props;
    cudaGetDevice(&device); cudaGetDeviceProperties(&device_props, device);
    
    // Determine grid and block sizes for kernel
    // NOTE : num_edges will never be more than the maz grid size of dim 0
    dim3 grid_size( num_edges, _snps / BLOCK_SIZE + 1, 1);
    dim3 block_size( 1, _snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, 1);
    
    // Invoke the search kernel 
    search_graph<<<grid_size, block_size, sizeof(size_t) * 2 * block_size.y>>>(_data_gpu, _graph, block_size.y);    
    CudaCheckError();
    
    // Sort the edges
    sort_edges(grid_size, block_size);
    
    // Create partitions
    map_to_partitions<<<grid_size, block_size>>>(_graph);
}

template <typename SubBlockType>
void Graph<SubBlockType, devices::gpu>::sort_edges(dim3 grid_size, dim3 block_size)
{
    const size_t    edges             = grid_size.x;
    const size_t    passes            = static_cast<size_t>(std::ceil(std::log2(static_cast<double>(edges)))); 
    size_t          out_in_block_size = 2;
    size_t          out_out_block_size;
    
    // Only need half the threads
    grid_size.x /= 2;
    
    for (size_t pass = 0; pass < passes; ++pass) {
        bitonic_out_in_sort<<<grid_size, block_size>>>(_graph, out_in_block_size, edges); 
        CudaCheckError();
        cudaThreadSynchronize();
        
        out_out_block_size = out_in_block_size / 2;
        for (size_t i = 0; i < pass; ++i) {
            bitonic_out_out_sort<<<grid_size, block_size>>>(_graph, out_out_block_size, edges);
            CudaCheckError();
            cudaThreadSynchronize();
            out_out_block_size /= 2;
        }
        out_in_block_size *= 2;
    }
    cudaThreadSynchronize();
#ifdef DEBUG 
    print_edges<<<1, 1>>>(_data_gpu, _graph);
    CudaCheckError();
#endif
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


