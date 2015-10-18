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
    bound_type*                 _selection_parameters;
public:    
    //-------------------------------------------------------------------------------------------------------
    /// @brief   Constructor 
    //-------------------------------------------------------------------------------------------------------
    CUDA_H
    TreeGpu(binary_vector& data  , read_info_container& read_info, snp_info_container         , 
            const size_t   snps  , const size_t         reads    , const size_t min_ubound    ,
            const size_t   device);

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Destuctor 
    //-------------------------------------------------------------------------------------------------------
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
        cudaFree(_selection_parameters);
    }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Solves the haplotype for the tree
    //-------------------------------------------------------------------------------------------------------
    void search_tree();
private:
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Allocates and copies memory to the device
    //-------------------------------------------------------------------------------------------------------
    void move_data();
    
};

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_GPU_H


