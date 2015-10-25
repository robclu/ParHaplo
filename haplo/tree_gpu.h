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
#define MAX_HEIGHT 21

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
    SubBlockType&               _sub_block;
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
    size_t                      _nih_cols;
    size_t                      _last_unaligned_read;
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
    Tree(SubBlockType& sub_block, const size_t device);

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
private:
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Puts the haplotype back into the sub_block 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H 
    void set_sub_block_haplotypes();
};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

template <typename SubBlockType>
Tree<SubBlockType, devices::gpu>::Tree(SubBlockType& sub_block, const size_t device)
: _sub_block(sub_block)                     ,
  _data(sub_block.data())                   , _read_info(sub_block.read_info())         , 
  _snp_info(sub_block.snp_info())           , _haplotype(sub_block.snp_info().size())   , 
  _alignments(sub_block.read_info().size())                                             , 
  _tree(sub_block.snp_info().size()         , sub_block.read_info().size())             , 
  _snps(sub_block.snp_info().size())        , _reads(sub_block.read_info().size())      ,
  _min_ubound(sub_block.size())             , _device(device)                           ,
  _nih_cols(sub_block.nih_columns())        , _last_unaligned_read(0)            
{
    cudaError_t error;
    
    // Copy the actual data to the device
    thrust::host_vector<small_type> data_gpu_format = _data.to_binary_vector();
    CudaSafeCall( cudaMalloc((void**)&_tree.data, sizeof(small_type) * _data.size()) );
    CudaSafeCall( cudaMemcpy(_tree.data, thrust::raw_pointer_cast(&data_gpu_format[0])   , 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice) );

    // Copy the read data to the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.read_info, sizeof(read_info_type) * _reads) );
    CudaSafeCall( cudaMemcpy(_tree.read_info, thrust::raw_pointer_cast(&_read_info[0]) ,
                    sizeof(read_info_type) * _reads, cudaMemcpyHostToDevice) );
    
    // Copy SNP data to the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.snp_info, sizeof(snp_info_type) * _snps) );
    CudaSafeCall( cudaMemcpy(_tree.snp_info, thrust::raw_pointer_cast(&_snp_info[0]) , 
                    sizeof(snp_info_type) * _snps, cudaMemcpyHostToDevice) );
  
    // Allocate the haplotype and alignment data on the device 
    CudaSafeCall( cudaMalloc((void**)&_tree.haplotype, sizeof(small_type) * 2 * _snps) );
    CudaSafeCall( cudaMemset(_tree.haplotype, 0, sizeof(small_type) * 2 * _snps)       );
    CudaSafeCall( cudaMalloc((void**)&_tree.alignments, sizeof(small_type) * _reads)   );
    CudaSafeCall( cudaMemset(_tree.alignments, 0, sizeof(small_type) *  _reads)        );
    
    // Create a vector for the search snps
    thrust::host_vector<size_t> snp_indices(_snps);
    thrust::sequence(snp_indices.begin(), snp_indices.end());
    
    // Copy snp search info to the device
    CudaSafeCall( cudaMalloc((void**)&_tree.search_snps, sizeof(size_t) * _snps) );
    CudaSafeCall( cudaMemcpy(_tree.search_snps, thrust::raw_pointer_cast(&snp_indices[0]), 
                    sizeof(size_t) * snp_indices.size(), cudaMemcpyHostToDevice)         );
    
    // Create a vector for the aligned rows
    thrust::host_vector<size_t> aligned_reads(_reads);
    thrust::sequence(aligned_reads.begin(), aligned_reads.end());
   
    // Copy aligned data to the device
    CudaSafeCall( cudaMalloc((void**)&_tree.aligned_reads, sizeof(size_t) * _reads) );
    CudaSafeCall( cudaMemcpy(_tree.aligned_reads, thrust::raw_pointer_cast(&aligned_reads[0]), 
                    sizeof(size_t) * aligned_reads.size(), cudaMemcpyHostToDevice)           );
    
    // Copy the array for the aligned offsets
    CudaSafeCall( cudaMalloc((void**)&_tree.alignment_offsets, sizeof(size_t) * _reads * 2) );
    CudaSafeCall( cudaMemset(_tree.alignment_offsets, 0, sizeof(size_t) * _reads * 2) );

    // Create data for bounds array
    CudaSafeCall( cudaMalloc((void**)&_snp_bounds, sizeof(BoundsGpu) * _snps) );

    // Allocate memory for the kernel arguments
    CudaSafeCall( cudaMalloc((void**)&_last_unaligned_idx, sizeof(size_t))  );
    CudaSafeCall( cudaMemcpy(_last_unaligned_idx, &_last_unaligned_read, 
                    sizeof(size_t), cudaMemcpyHostToDevice            ));
    
    // Offset in the alignment array
    CudaSafeCall( cudaMalloc((void**)&_alignment_offset, sizeof(size_t))  );
    CudaSafeCall( cudaMemset(_alignment_offset, 0, sizeof(size_t))        );

    size_t free_mem, total_mem;                 // Device memory
    
    // --------------------------------------- READ VALUES --------------------------------------------------
    
    // If the tree is small, we can allocate the "correct" (worst case) amount of memory
    if (_snps <= BIG_INPUT) { 
        CudaSafeCall( cudaMalloc((void**)&_tree.read_values, sizeof(uint8_t) * pow(2, _snps)) );
        CudaSafeCall( cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * pow(2, _snps)) );
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
    std::cout << "NIH COLUMNS : " << _nih_cols << "\n";
    
    size_t unsearched_snps   = _snps - 1;            // SNPs still to search 
    size_t prev_level_start  = 0;                    // Start index in array for previous level
    size_t this_level_start  = 1;                    // Start index in array for this level
    size_t nodes_in_level    = 2;                    // Number of nodes in the level
    size_t last_searched_snp = 0;                    // Index of the last searched snp
    int    device;
    
    cudaDeviceProp device_props;                    // Properties of the device
   
    // ---------------------------------- DEVICE PROPERTIES -------------------------------------------------
    
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&device_props, device); // Get the device properties
    
    // --------------------------------------- ROOT NODE ----------------------------------------------------
    
    // The first level of the tree (with the root node) is a little different
    // because it needs to do the setup, so invoke that kernel first
    map_root_node<<<1, 1>>>(_tree            , _snp_bounds, last_searched_snp, _last_unaligned_idx, 
                            _alignment_offset, _min_ubound, static_cast<size_t>(device)          ); 
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

    _nih_cols -= 35;
    size_t terminate = 0;
    while (last_searched_snp < _snps - _nih_cols) {
        
        // We need to call the grid manager here 
        dim3 grid_size( (nodes_in_level / BLOCK_SIZE + 1) % device_props.maxGridSize[0]    , 
                        (nodes_in_level / BLOCK_SIZE + 1) / device_props.maxGridSize[0] + 1, 
                        1);
        dim3 block_size( nodes_in_level < BLOCK_SIZE ? (nodes_in_level / WARP_SIZE + 1) * WARP_SIZE : BLOCK_SIZE,
                         grid_size.y > 1 ? BLOCK_SIZE : 1,
                         1);
        
        // Perform a "mapping" step, which maps the bounds onto the nodes in the level
        map_level<<<grid_size, block_size>>>(_tree, _snp_bounds, last_searched_snp, _last_unaligned_idx, 
                                             _alignment_offset , prev_level_start , this_level_start   , 
                                             nodes_in_level    );
        CudaCheckError();

        // "Reduce" the search space to eliminate the bad nodes
        reduce_level<<<grid_size, block_size>>>(_tree, this_level_start, nodes_in_level);
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
        
        if (++terminate < MAX_HEIGHT)
            nodes_in_level   *= 2;
    }
   
    // Get the number of unaligned reads 
    CudaSafeCall( cudaMemcpy(&_last_unaligned_read, _last_unaligned_idx, sizeof(size_t), cudaMemcpyDeviceToHost) );    
    CudaCheckError();
    
    // Create a block for each of the unaligned snps
    dim3 grid_size(_reads - _last_unaligned_read, 1, 1);
    dim3 block_size(_nih_cols > BLOCK_SIZE ? BLOCK_SIZE : _nih_cols, _nih_cols / BLOCK_SIZE + 1, 1);
    
    // Align all the reads which were not aligned
    align_unaligned_reads<<<grid_size, block_size, _nih_cols * 2 * sizeof(size_t)>>>(_tree, _last_unaligned_idx,
                                                                                     last_searched_snp);
    
//    is_sorted<<<1, 1>>>(_tree, prev_level_start, nodes_in_level / 2);
    // Find the haplotype for the ih columns
    find_temp_haplotype<<<1, 1>>>(_tree, prev_level_start);
   
    // Create a block for each of the unsolved nih columns
    grid_size  = dim3(_snps - last_searched_snp - 1, 1, 1);
    block_size = dim3(_reads > BLOCK_SIZE ? BLOCK_SIZE : _reads, _reads / BLOCK_SIZE + 1, 1);
    
    // Lastly, find the final haplotype which includes the nih columns
    solve_nih_columns<<<grid_size, block_size, _reads * 3 * sizeof(size_t)>>>(_tree, last_searched_snp);
    cudaDeviceSynchronize();
    
    // Move the haplotypes to the sub-block 
    set_sub_block_haplotypes();
}

template <typename SubBlockType>
void Tree<SubBlockType, devices::gpu>::set_sub_block_haplotypes()
{
   // Create a host vector to get the haplotype data into
    thrust::host_vector<uint8_t> haplotypes(_snps * 2);
    
    // Get the data 
    CudaSafeCall( cudaMemcpy(thrust::raw_pointer_cast(&haplotypes[0]), _tree.haplotype,
                        sizeof(uint8_t) * _snps * 2, cudaMemcpyDeviceToHost           ));
    
    // Go trough the haplotypes and set the values in the sub_block
    for (size_t i = 0; i < _snps; ++i) {
        _sub_block._haplo_one.set(i, haplotypes[i]);
        _sub_block._haplo_two.set(i, haplotypes[i + _snps]);
    }
}


}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


