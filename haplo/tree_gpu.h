// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_GPU_H
#define PARHAPLO_TREE_GPU_H

#include "cuda_error.h"
#include "devices.hpp"
#include "grid_config.h"
#include "tree.hpp"
#include "tree_kernels_gpu.cu"
#include "thrust/sequence.h"

#define BIG_INPUT  24
#define ALIGN_MEM  0.3f
#define NODES_MEM  0.9f

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
    size_t*                     _last_unaligned_idx;
    size_t*                     _alignment_offset;
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
        cudaFree(_tree.data)                ;
        cudaFree(_tree.read_info)           ;
        cudaFree(_tree.snp_info)            ;
        cudaFree(_tree.haplotype)           ;
        cudaFree(_tree.alignments)          ;
        cudaFree(_tree.search_snps)         ;
        cudaFree(_tree.aligned_reads)       ;
        cudaFree(_tree.alignment_offsets)   ;
        cudaFree(_snp_bounds)               ;
        cudaFree(_tree.read_values)         ;
        cudaFree(_tree.nodes)               ;
        cudaFree(_last_unaligned_idx)       ;
        cudaFree(_alignment_offset)         ;
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
    
    // Copy the actual data to the device
    thrust::host_vector<small_type> data_gpu_format = data.to_binary_vector();
    CudaSafeCall( cudaMalloc((void**)&_tree.data, sizeof(small_type) * data.size()) );
    CudaSafeCall( cudaMemcpy(_tree.data, thrust::raw_pointer_cast(&data_gpu_format[0])   , 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice) );

    // Copy the read data to the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.read_info, sizeof(read_info_type) * read_info.size()) );
    CudaSafeCall( cudaMemcpy(_tree.read_info, thrust::raw_pointer_cast(&_read_info[0]) ,
                    sizeof(read_info_type) * read_info.size(), cudaMemcpyHostToDevice) );
    
    // Copy SNP data to the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.snp_info, sizeof(snp_info_type) * snp_info.size()) );
    CudaSafeCall( cudaMemcpy(_tree.snp_info, thrust::raw_pointer_cast(&_snp_info[0]) , 
                    sizeof(snp_info_type) * snp_info.size(), cudaMemcpyHostToDevice) );
  
    // Allocate the haplotype and alignment data on the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.haplotype, sizeof(small_type) * snps) );
    CudaSafeCall( cudaMalloc((void**)&_tree.alignments, sizeof(small_type) * read_info.size()) );;
    
    // Create a vector for the search snps
    thrust::host_vector<size_t> snp_indices(snps);
    thrust::sequence(snp_indices.begin(), snp_indices.end());
    
    // Copy snp search info to the device
    CudaSafeCall( cudaMalloc((void**)&_tree.search_snps, sizeof(size_t) * snps) );
    CudaSafeCall( cudaMemcpy(_tree.search_snps, thrust::raw_pointer_cast(&snp_indices[0]), 
                    sizeof(size_t) * snp_indices.size(), cudaMemcpyHostToDevice)         );
    
    // Create a vector for the aligned rows
    thrust::host_vector<size_t> aligned_reads(reads);
    thrust::sequence(aligned_reads.begin(), aligned_reads.end());
   
    // Copy aligned data to the device
    CudaSafeCall( cudaMalloc((void**)&_tree.aligned_reads, sizeof(size_t) * reads) );
    CudaSafeCall( cudaMemcpy(_tree.aligned_reads, thrust::raw_pointer_cast(&aligned_reads[0]), 
                    sizeof(size_t) * aligned_reads.size(), cudaMemcpyHostToDevice)           );
    
    // Copy the array for the aligned offsets
    CudaSafeCall( cudaMalloc((void**)&_tree.alignment_offsets, sizeof(size_t) * reads * 2) );
    CudaSafeCall( cudaMemset(_tree.alignment_offsets, 0, sizeof(size_t) * reads * 2) );

    // Create data for bounds array
    CudaSafeCall( cudaMalloc((void**)&_snp_bounds, sizeof(BoundsGpu) * snps) );

    // Allocate memory for the kernel arguments
    CudaSafeCall( cudaMalloc((void**)&_last_unaligned_idx, sizeof(size_t))  );
    CudaSafeCall( cudaMemset(_last_unaligned_idx, 0, sizeof(size_t))        );
    
    // Offset in the alignment array
    CudaSafeCall( cudaMalloc((void**)&_alignment_offset, sizeof(size_t))  );
    CudaSafeCall( cudaMemset(_alignment_offset, 0, sizeof(size_t))        );

    size_t free_mem, total_mem;                 // Device memory
    
    // --------------------------------------- READ VALUES --------------------------------------------------
    
    // If the tree is small, we can allocate the "correct" (worst case) amount of memory
    if (snps <= BIG_INPUT) { 
        CudaSafeCall( cudaMalloc((void**)&_tree.read_values, sizeof(uint8_t) * pow(2, snps)) );
        CudaSafeCall( cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * pow(2, snps)) );
    } else {
        cudaMemGetInfo(&free_mem, &total_mem);
        // Allocate space for the alignments for the nodes
        CudaSafeCall( cudaMalloc((void**)&_tree.read_values, ALIGN_MEM * free_mem) );
        cudaMemGetInfo(&free_mem, &total_mem);
        // Allocate space for the nodes themselves
        size_t num_nodes = NODES_MEM * free_mem / sizeof(TreeNode);
        CudaSafeCall( cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * num_nodes) );
#ifdef DEBUG
        cudaMemGetInfo(&free_mem, &total_mem);
        std::cout << "Free  Memory : " << free_mem  << "\n";
        std::cout << "Total Memory : " << total_mem << "\n";
        std::cout << "Total Nodes  : " << num_nodes << "\n";
#endif
    }
    
    // ---------------------------------------- HEAP LIMIT ---------------------------------------------------
    
    const size_t heap_mb = 64;                         // Excessive, but just incase
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
    map_root_node<<<1, 1>>>(_tree            , _snp_bounds, last_searched_snp, _last_unaligned_idx, 
                            _alignment_offset, _min_ubound, _device          ); 
    CudaCheckError();

    // Now we can do the mapping of the unsearched nodes
    map_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, last_searched_snp, 0);
    CudaCheckError();

    // And the reduction
    reduce_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, last_searched_snp, unsearched_snps); 
    CudaCheckError();

    ++last_searched_snp;    
    --unsearched_snps;
    
    // ----------------------------------------- OTHER NODES ------------------------------------------------

    size_t terminate = 0;
    while (last_searched_snp < _snps && terminate++ < 19) {
        
        // We need to call the grid manager here 
        dim3 grid_size(nodes_in_level / 1024 + 1, 1, 1);
        const size_t threads = nodes_in_level < 1024 
                             ? (nodes_in_level / 32 + 1) * 32 : 1024;
        
        // Perform a "mapping" step, which maps the bounds onto the nodes in the level
        map_level<<<grid_size, threads>>>(_tree, _snp_bounds, last_searched_snp, _last_unaligned_idx, 
                                          _alignment_offset , prev_level_start , this_level_start   , 
                                          nodes_in_level    );
        CudaCheckError();

        // "Reduce" the search space to eliminate the bad nodes
        reduce_level<<<grid_size, threads>>>(_tree, this_level_start, nodes_in_level);
        CudaCheckError();

        // Map unsearched snps
        map_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, last_searched_snp, this_level_start);
        CudaCheckError();
        
        // Reduce unsearched snps        
        reduce_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, last_searched_snp, unsearched_snps);
        CudaCheckError();

        --unsearched_snps;
        ++last_searched_snp;
        
        // Update variables for next iteration 
        prev_level_start  = this_level_start;
        this_level_start += nodes_in_level;
        nodes_in_level   *= 2;
    }
   
    // Haplotype found
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


