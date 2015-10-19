// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_GPU_H
#define PARHAPLO_TREE_GPU_H

#include "cuda_error.h"
#include "devices.hpp"
#include "tree.hpp"
#include "tree_kernels_gpu.cu"

namespace haplo {

// Dummy block class since this doesn't actually need the sub_block as a template
class VoidSubBlock;

// Specialization for gpu
template <typename SubBlockType>    
class Tree<SubBlockType, devices::gpu> {
public:
    //-------------------------------------------------------------------------------------------------------
    using binary_vector                 = BinaryVector<2>;
    using tree_type                     = internal::Tree;
    using read_info_type                = typename tree_type::read_info_type;
    using read_info_container           = thrust::host_vector<read_info_type>;
    using snp_info_type                 = typename tree_type::snp_info_type;
    using snp_info_container            = thrust::host_vector<snp_info_type>;
    using small_type                    = typename tree_type::small_type;
    using small_container               = thrust::host_vector<small_type>;
    using bound_type                    = BoundsGpu;
    //-------------------------------------------------------------------------------------------------------
private:
    // -------------------------------------------- HOST ----------------------------------------------------
    binary_vector&              _data;
    read_info_container&        _read_info;
    snp_info_container          _snp_info;
    small_container             _haplotype;
    small_container             _alignments;
    size_t                      _snps;
    size_t                      _reads;
    size_t                      _min_ubound;
    size_t                      _device;
    // ------------------------------------------ DEVICE ----------------------------------------------------
    tree_type                   _tree; 
    bound_type*                 _snp_bounds;
    size_t*                     _node_index;
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    Tree(binary_vector& data  , read_info_container& read_info, snp_info_container         , 
         const size_t   snps  , const size_t         reads    , const size_t min_ubound    ,
         const size_t   device);

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Destuctor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    ~Tree() 
    {
        cudaFree(_tree.data);
        cudaFree(_tree.read_info);
        cudaFree(_tree.snp_info);
        cudaFree(_tree.nodes);
        cudaFree(_tree.haplotype);
        cudaFree(_tree.alignments);
        cudaFree(_tree.search_snps);
        cudaFree(_tree.aligned_reads);
        cudaFree(_snp_bounds);
        cudaFree(_node_index);
    }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Solves the haplotype for the tree
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    void search();
};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

template <typename SubBlockType>
Tree<SubBlockType, devices::gpu>::Tree(
                 binary_vector& data    , read_info_container& read_info, snp_info_container snp_info   ,   
                 const size_t   snps    , const size_t         reads    , const size_t       min_ubound ,
                 const size_t   device  )
: _data(data)            , _read_info(read_info) , _snp_info(snp_info)    , _haplotype(snps)  , 
  _alignments(reads)     , _tree(snps, reads)    , _snps(snps)            , _reads(reads)     ,
  _min_ubound(min_ubound), _device(device)
{
    cudaError_t error;
    
    // Copy data to device
    thrust::host_vector<small_type> data_gpu_format = data.to_binary_vector();
    CudaSafeCall( cudaMalloc((void**)&_tree.data, sizeof(small_type) * data.size()) );
    CudaSafeCall( cudaMemcpy(_tree.data, thrust::raw_pointer_cast(&data_gpu_format[0])   , 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice) );

    // Copy the read data to the device 
    cudaMalloc((void**)&_tree.read_info, sizeof(read_info_type) * read_info.size());
    cudaMemcpy(_tree.read_info, thrust::raw_pointer_cast(&_read_info[0]),
                    sizeof(read_info_type) * read_info.size(), cudaMemcpyHostToDevice); 
    
    // Copy SNP data to the device 
    cudaMalloc((void**)&_tree.snp_info, sizeof(snp_info_type) * snp_info.size());
    cudaMemcpy(_tree.snp_info, thrust::raw_pointer_cast(&_snp_info[0]), 
                    sizeof(snp_info_type) * snp_info.size(), cudaMemcpyHostToDevice); 
  
    // Allocate the haplotype and alignment data on the device 
    cudaMalloc((void**)&_tree.haplotype, sizeof(small_type) * snps);
    cudaMalloc((void**)&_tree.alignments, sizeof(small_type) * read_info.size());
    
    // Create a vector for the search snps
    thrust::host_vector<size_t> snp_indices;
    for (size_t i = 0; i < snps; ++i) snp_indices.push_back(i);
    
    // Copy snp search info to the device
    cudaMalloc((void**)&_tree.search_snps, sizeof(size_t) * snps);
    cudaMemcpy(_tree.search_snps, thrust::raw_pointer_cast(&snp_indices[0]), 
                sizeof(size_t) * snp_indices.size(), cudaMemcpyHostToDevice);
    
    // Create a vector for the aligned rows
    thrust::host_vector<size_t> aligned;
    for (size_t i = 0; i < reads; ++i) aligned.push_back(i);
   
    // Copy aligned data to the device
    cudaMalloc((void**)&_tree.aligned_reads, sizeof(size_t) * reads);
    cudaMemcpy(_tree.aligned_reads, thrust::raw_pointer_cast(&aligned[0]), 
                sizeof(size_t) * aligned.size(), cudaMemcpyHostToDevice);
    
    // Create data for bounds array
    cudaMalloc((void**)&_snp_bounds, sizeof(BoundsGpu) * snps);

    // Allocate memory for the kernel arguments
    CudaSafeCall( cudaMalloc((void**)&_node_index, sizeof(size_t))  );
    CudaSafeCall( cudaMemset(_node_index, 0, sizeof(size_t))        );

    // ------------------------------------ NODE ALLOCATION -------------------------------------------------
    
    // Allocate space on the device for the tree nodes
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
   
    // DEBUG 
    std::cout << "Free Memory Before : " << free_mem << "\n";
    
    // Determine the number of nodes to "safely allocate"
    size_t num_nodes = 0.7 * free_mem / sizeof(TreeNode);
    
    // Check that we can allocate the memory
    error = cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * num_nodes);
    if (error != cudaSuccess) std::cout << "Error allocating memory\n";
    
    // DEBUG
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "Free Memory After: " << free_mem  << "\n";
    std::cout << "Num Nodes        : " << num_nodes << "\n";
    
    // ------------------------------------- KERNEL ARGUMENTS -----------------------------------------------
    

    // ---------------------------------------- HEAP LIMIT ---------------------------------------------------
    
    const size_t heap_mb = 128;                         // Excessive, but just incase
    cudaThreadSetLimit(cudaLimitMallocHeapSize, heap_mb * 1024 * 1024);
}

template <typename SubBlockType>
void Tree<SubBlockType, devices::gpu>::search()
{
    size_t unsearched_snps   = _snps - 1;            // SNPs still to search 
    size_t prev_level_start  = 0;                    // Start index in array for previous level
    size_t this_level_start  = 1;                    // Start index in array for this level
    size_t nodes_in_level    = 2;                    // Number of nodes in the level
    size_t last_searched_snp = 0;                    // Index of the last searched snp
   
    // --------------------------------------- ROOT NODE ----------------------------------------------------
    
    // The first level of the tree (with the root node) is a little different
    // because it needs to do the setup, so invoke that kernel first
    map_root_node<<<1, 1>>>(_tree, _snp_bounds, _min_ubound, _device); 
    CudaCheckError();
    
    // Now we can do the mapping of the unsearched nodes
    map_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, _node_index);
    CudaCheckError();
    
/*    
    // And the reduction
    reduce_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, unsearched_snps); 
    if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Unsearched SNP Reduce!\n");
    cudaDeviceSynchronize();
    
    --unsearched_snps;
    
    // ----------------------------------------- OTHER NODES ------------------------------------------------

    size_t terminate = 0;
    while (last_searched_snp < _snps && terminate < 4) {
        // Perform a "mapping" step, which maps the bounds onto the nodes in the level
        map_level<<<1, nodes_in_level>>>(_tree, _snp_bounds, prev_level_start, this_level_start);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Level Map!\n");
        cudaDeviceSynchronize();
       
        // Need to add in level reduction here to reduce search space 
        
        // Map unsearched snps
        map_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, this_level_start);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Unsearched SNP Map!\n");
        cudaDeviceSynchronize();

        // Reduce unsearched snps        
        reduce_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, unsearched_snps);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Unsearched SNP Reduce!\n");
        cudaDeviceSynchronize();
        
        --unsearched_snps;
        
        // Update variables for next iteration 
        prev_level_start  = this_level_start;
        this_level_start += nodes_in_level;
        nodes_in_level   *= 2;
    }
  */  
    // Haplotype found
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


