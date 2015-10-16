// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class -- gpu implementattion
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_GPU_H
#define PARHAPLO_TREE_GPU_H

#include "tree_kernels_gpu.cu"

namespace haplo {
    
class TreeGpu {
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
    //-------------------------------------------------------------------------------------------------------
private:
    // -------------------------------------------- HOST ----------------------------------------------------
    binary_vector&              _data;
    read_info_container&        _read_info;
    snp_info_container          _snp_info;
    small_container             _haplotype;
    small_container             _alignments;
    // ------------------------------------------ DEVICE ----------------------------------------------------
    tree_type                   _tree; 
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    TreeGpu(binary_vector& data, read_info_container& read_info, snp_info_container , 
            const size_t   snps, const size_t         reads                         );
    
    CUDA_H
    ~TreeGpu() 
    {
        cudaFree(_tree.data);
        cudaFree(_tree.read_info);
        cudaFree(_tree.snp_info);
        cudaFree(_tree.nodes);
        cudaFree(_tree.haplotype);
        cudaFree(_tree.alignments);
        cudaFree(_tree.search_snps);
        cudaFree(_tree.aligned_reads);
    }
};


// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

TreeGpu::TreeGpu(binary_vector& data, read_info_container& read_info, snp_info_container snp_info, 
                 const size_t snps, const size_t reads)
: _data(data)       , _read_info(read_info) , _snp_info(snp_info),
  _haplotype(snps)  , _alignments(reads)    , _tree(snps, reads)
{
    cudaError_t error;
    
    // Copy data to device
    thrust::host_vector<small_type> data_gpu_format = data.to_binary_vector();
    cudaMalloc((void**)&_tree.data, sizeof(small_type) * data.size());
    cudaMemcpy(_tree.data, thrust::raw_pointer_cast(&data_gpu_format[0]), 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice); 

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


    // ------------------------------------ NODE ALLOCATION -------------------------------------------------
    
    // Allocate space on the device for the tree nodes
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
   
    // DEBUG 
    std::cout << "Free Memory Before : " << free_mem << "\n";
    
    // Determine the number of nodes to "safely allocate"
    size_t num_nodes                 = 0.1 * free_mem / sizeof(TreeNode);
    
    // Check that we can allocate the memory
    error = cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * num_nodes);
    if (error != cudaSuccess) std::cout << "Error allocating memory\n";
    
    // DEBUG
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "Free Memory After: " << free_mem << "\n";

    // --------------------------------------- TREE SEARCH ---------------------------------------------------

    // Invoke the kernel with just a single thread -- kernel spawns more threads
    search_tree<<<1,1>>>(_tree);
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


