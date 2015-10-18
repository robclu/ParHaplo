// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for internal tree sruct to clean the GPU interface
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_INTERNAL_H
#define PARHAPLO_TREE_INTERNAL_H

#include "read_info_gpu.cuh"
#include "snp_info_gpu.h"
#include "tree_node.h"
#include "thrust/host_vector.h"

namespace haplo     {
namespace internal  {
 
struct Tree {    
public:
    //-------------------------------------------------------------------------------------------------------
    using read_info_type                = ReadInfo;
    using snp_info_type                 = SnpInfoGpu;
    using small_type                    = uint8_t;
    using node_type                     = TreeNode;
    using small_container               = thrust::host_vector<small_type>;
    //------------------------------------------------------------------------------------------------------- 
   
    //------------------------------------------------------------------------------------------------------- 
    /// @brief  Constructor to initialize avriables
    //------------------------------------------------------------------------------------------------------- 
    Tree(const size_t num_snps, const size_t num_reads)
    : snps(num_snps), reads(num_reads), last_searched_snp(0), last_unaligned_idx(0) {}

    size_t          snps;
    size_t          reads;
    size_t          last_searched_snp;
    size_t          last_unaligned_idx;
    size_t*         search_snps;
    size_t*         aligned_reads;
    small_type*     data;
    small_type*     haplotype;
    small_type*     alignments;
    read_info_type* read_info;
    snp_info_type*  snp_info;
    node_type*      nodes; 
    
    //------------------------------------------------------------------------------------------------------- 
    /// @brief      Gets a pointer to one of the nodes
    //------------------------------------------------------------------------------------------------------- 
    CUDA_HD
    const TreeNode* node_ptr(const size_t i) const { return (const TreeNode*)(&nodes[i]); }
    
    //------------------------------------------------------------------------------------------------------- 
    /// @brief      Gets a pointer to one of the nodes
    //------------------------------------------------------------------------------------------------------- 
    CUDA_HD
    TreeNode* node_ptr(const size_t i) { return &nodes[i]; }
};

}
}
#endif          // PARAHAPLO_TREE_INTERNAL_H
