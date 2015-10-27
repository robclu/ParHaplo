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
#define ITERS            10

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
    size_t                      _mec_score;
    
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
        cudaFree(_data_gpu.data       );
        cudaFree(_data_gpu.read_info  );
        cudaFree(_data_gpu.snp_info   );
        cudaFree(_graph.edges         );
        cudaFree(_graph.set_one       );
        cudaFree(_graph.set_two       );
        cudaFree(_graph.set_one       );
        cudaFree(_graph.set_two       );
        cudaFree(_graph.haplo_one     );
        cudaFree(_graph.haplo_two     );
        cudaFree(_graph.haplo_one_temp);
        cudaFree(_graph.haplo_two_temp);
        cudaFree(_graph.fragments     );
        cudaFree(_graph.mec_score     );
    }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Solves the graph for the haplotypes
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void search();
private:
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Sorts the edges
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void sort_edges(dim3 grid_size, dim3 block_size);

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Sorts the fragments
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void sort_fragments(dim3 grid_size, dim3 block_size);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Refines the solution until the best score is reached
    // ------------------------------------------------------------------------------------------------------
    CUDA_H 
    size_t refine_solution(const cudaStream_t* streams, const size_t mem_size);
    
    CUDA_H 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Moves the result of the haplotype to the sub block 
    // ------------------------------------------------------------------------------------------------------
    void set_sub_block_haplotypes();

};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

template <typename SubBlockType>
Graph<SubBlockType, devices::gpu>::Graph(SubBlockType& sub_block, const size_t device)
: _sub_block(sub_block)                     , _data_cpu(sub_block.data())               , 
  _read_info(sub_block.read_info())         , _snp_info(sub_block.snp_info())           , 
  _haplotype(sub_block.snp_info().size())   , _alignments(sub_block.read_info().size()) , 
  _snps(sub_block.snp_info().size())        , _reads(sub_block.read_info().size())      ,
  _device(device)                           , _nih_cols(sub_block.nih_columns())        ,
  _mec_score(INT_MAX)                       ,                   
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
    
    // Allocate space for the partitions
    CudaSafeCall( cudaMalloc((void**)&_graph.set_one, sizeof(size_t) * _reads) );
    CudaSafeCall( cudaMalloc((void**)&_graph.set_two, sizeof(size_t) * _reads) );
    CudaSafeCall( cudaMemset(_graph.set_one, INT_MAX, sizeof(size_t) * _reads) );
    CudaSafeCall( cudaMemset(_graph.set_two, INT_MAX, sizeof(size_t) * _reads) );
    
    // Allocate space for the haplotypes
    CudaSafeCall( cudaMalloc((void**)&_graph.haplo_one, sizeof(size_t) * _snps)      );
    CudaSafeCall( cudaMalloc((void**)&_graph.haplo_two, sizeof(size_t) * _snps)      );
    CudaSafeCall( cudaMalloc((void**)&_graph.haplo_one_temp, sizeof(size_t) * _snps) );
    CudaSafeCall( cudaMalloc((void**)&_graph.haplo_two_temp, sizeof(size_t) * _snps) );
    
    // Allocate space for each of the values which contribute to the switch error rate (MEC score)
    CudaSafeCall( cudaMalloc((void**)&_graph.set_one_counts, sizeof(size_t) * _snps * 2) );
    CudaSafeCall( cudaMalloc((void**)&_graph.set_two_counts, sizeof(size_t) * _snps * 2) );
    CudaSafeCall( cudaMalloc((void**)&_graph.fragments, sizeof(Fragment) * _reads)       );
    
    // Allocate device memory for the mec score
    CudaSafeCall( cudaMalloc((void**)&_graph.mec_score, sizeof(size_t)) );
    CudaSafeCall( cudaMemset(_graph.mec_score, INT_MAX, sizeof(size_t)) );
 
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
    map_to_partitions<<<grid_size, block_size, sizeof(uint8_t) * _data_gpu.reads * 2 >>>(_data_gpu, _graph);
    CudaCheckError();
    cudaDeviceSynchronize();
    
    // Create streams 
    cudaStream_t streams[2];
    cudaStreamCreate(&streams[0]); cudaStreamCreate(&streams[1]);
    
    // Resize the grid 
    grid_size = dim3(_snps, _reads / BLOCK_SIZE + 1, 1);
    block_size = dim3(1, _reads > BLOCK_SIZE ? BLOCK_SIZE : _reads, 1);

    const size_t mem_size = sizeof(size_t) * _reads * 2;    
    
    // Determine the starting score for the haplotype
    determine_switch_error<1><<<grid_size, block_size, mem_size, streams[0]>>>(_data_gpu, _graph);
    CudaCheckError();
    
    determine_switch_error<2><<<grid_size, block_size, mem_size, streams[1]>>>(_data_gpu, _graph); 
    CudaCheckError();
    cudaDeviceSynchronize();
   
    // Partition the unpartitioned reads
    grid_size = dim3(_reads, _snps / BLOCK_SIZE + 1, 1);
    block_size = dim3(1, _snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, 1);
    
    add_unpartitioned<<<grid_size, block_size, sizeof(size_t) * 2 * _snps>>>(_data_gpu, _graph);
    CudaCheckError();
   
    // Maps the score of each of the fragments based on the current haplotypes
    map_mec_score<<<grid_size, block_size, sizeof(size_t) * 2 * _snps>>>(_data_gpu, _graph);
    CudaCheckError();
    
    // Reduces the fragment MEC scores to get the overall MEC score
    block_size = dim3(_snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, _snps / BLOCK_SIZE + 1, 1);
    reduce_mec_score<<<1, block_size, sizeof(size_t) * _reads>>>(_data_gpu, _graph);
    CudaCheckError();

    // Refine the solution
    size_t prev_mec_score, terminate = 0;
    do {
        prev_mec_score = refine_solution(&streams[0], mem_size);
        if (prev_mec_score == _mec_score) ++terminate;
    } while (prev_mec_score >= _mec_score && terminate < ITERS);
    
    // Print the haplotypes 
    print_haplotypes<<<1, 1>>>(_data_gpu, _graph);
    CudaCheckError();
    
    // Put the haplotypes back into the sub_block 
    set_sub_block_haplotypes();
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
}

template <typename SubBlockType>
void Graph<SubBlockType, devices::gpu>::sort_fragments(dim3 grid_size, dim3 block_size)
{
    const size_t    fragments         = grid_size.x;
    const size_t    passes            = static_cast<size_t>(std::ceil(std::log2(static_cast<double>(fragments)))); 
    size_t          out_in_block_size = 2;
    size_t          out_out_block_size;
    
    // Only need half the threads
    grid_size.x /= 2;
    
    for (size_t pass = 0; pass < passes; ++pass) {
        bitonic_out_in_sort_frag<<<grid_size, block_size>>>(_graph, out_in_block_size, fragments); 
        CudaCheckError();
        cudaThreadSynchronize();
        
        out_out_block_size = out_in_block_size / 2;
        for (size_t i = 0; i < pass; ++i) {
            bitonic_out_out_sort_frag<<<grid_size, block_size>>>(_graph, out_out_block_size, fragments);
            CudaCheckError();
            cudaThreadSynchronize();
            out_out_block_size /= 2;
        }
        out_in_block_size *= 2;
    }
    cudaThreadSynchronize();
}

template <typename SubBlockType>
size_t Graph<SubBlockType, devices::gpu>::refine_solution(const cudaStream_t* streams  , 
                                                          const size_t        mem_size )
{
    size_t mec_score_before = _mec_score;
    
    dim3 grid_size(_reads, _snps / BLOCK_SIZE + 1, 1);
    dim3 block_size(_snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, _snps / BLOCK_SIZE + 1, 1);
    
    // Determine switch error rate of the current solution -- make a loop
    determine_switch_error<1><<<grid_size, block_size, mem_size, streams[0]>>>(_data_gpu, _graph);
    CudaCheckError();
    
    determine_switch_error<2><<<grid_size, block_size, mem_size, streams[1]>>>(_data_gpu, _graph); 
    CudaCheckError();
    cudaDeviceSynchronize();
/*
    // Check if flipping NIH column bits will give a better solution
    grid_size  = dim3(_snps, _reads / BLOCK_SIZE + 1, 1);
    block_size = dim3(1, _reads > BLOCK_SIZE ? BLOCK_SIZE : _reads, 1);
    
    evaluate_nih_columns<1><<<grid_size, block_size, mem_size, streams[0]>>>(_data_gpu, _graph);
    CudaCheckError();
    evaluate_nih_columns<2><<<grid_size, block_size, mem_size, streams[1]>>>(_data_gpu, _graph);
    CudaCheckError();
    cudaDeviceSynchronize();
*/
    grid_size = dim3(_reads, _snps / BLOCK_SIZE + 1, 1);
    block_size = dim3(1, _snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, 1);
    
    // Maps the score of each of the fragments based on the current haplotypes
    map_mec_score<<<grid_size, block_size, sizeof(size_t) * 2 * _snps>>>(_data_gpu, _graph);
    CudaCheckError();
    
    // Reduces the fragment MEC scores to get the overall MEC score
    block_size = dim3(_snps > BLOCK_SIZE ? BLOCK_SIZE : _snps, _snps / BLOCK_SIZE + 1, 1);
    reduce_mec_score<<<1, block_size, sizeof(size_t) * _reads>>>(_data_gpu, _graph);
    CudaCheckError(); 
    cudaDeviceSynchronize();
    
    // Copy the MEC score back 
    CudaSafeCall( cudaMemcpy(&_mec_score, _graph.mec_score, sizeof(size_t), cudaMemcpyDeviceToHost) );
    std::cout << "MEC SCORE : " << _mec_score << "\n";
    
    sort_fragments(dim3(_reads,1 ,1), dim3(1, 1, 1));
    CudaCheckError();
    swap_fragment_set<<<1, 1>>>(_data_gpu, _graph);
    CudaCheckError();
    
    return mec_score_before;
}

template <typename SubBlockType>
void Graph<SubBlockType, devices::gpu>::set_sub_block_haplotypes()
{
    // Create the host vectors 
    thrust::host_vector<uint8_t> haplo_one(_snps), haplo_two(_snps);
    
    // Get the results from the GPU
    CudaSafeCall( cudaMemcpy(thrust::raw_pointer_cast(&haplo_one[0]), _graph.haplo_one, 
                    sizeof(uint8_t) * _snps, cudaMemcpyDeviceToHost                   ) );
    CudaSafeCall( cudaMemcpy(thrust::raw_pointer_cast(&haplo_two[0]), _graph.haplo_two, 
                    sizeof(uint8_t) * _snps, cudaMemcpyDeviceToHost                   ) );
    
    for (size_t i = 0; i < _snps; ++i) {
        _sub_block._haplo_one.set(i, haplo_one[i]);
        _sub_block._haplo_two.set(i, haplo_two[i]);
    }
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


