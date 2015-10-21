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
    //------------------------------------- ALIAS's ---------------------------------------------------------
    using read_info_type                = ReadInfo;
    using snp_info_type                 = SnpInfoGpu;
    using small_type                    = uint8_t;
    using node_type                     = TreeNode;
    using small_container               = thrust::host_vector<small_type>;
    //------------------------------------------------------------------------------------------------------- 

    size_t          snps;               // Number of snps to serach
    size_t          reads;              // Number of reads for the sub-block
    size_t*         search_snps;        // Snps (indices) which have been (and still need to be) searched
    size_t*         aligned_reads;      // Reads (indices) which have been (and still need to be) searched
    size_t*         alignment_offsets;  // The offset of the start of each set of alignments for a level 
    small_type*     read_values;        // The values of the alignment of the reads for each node
    small_type*     data;               // The actual data 
    small_type*     haplotype;          // The final haplotype
    small_type*     alignments;         // The final alignments
    read_info_type* read_info;          // The information for each read for fast access
    snp_info_type*  snp_info;           // The information for each snp for fast access
    node_type*      nodes;              // The nodes which make up the tree

    //------------------------------------------------------------------------------------------------------- 
    /// @brief  Constructor to initialize avriables
    //------------------------------------------------------------------------------------------------------- 
    Tree(const size_t num_snps, const size_t num_reads)
    : snps(num_snps), reads(num_reads) {}

    //------------------------------------------------------------------------------------------------------- 
    /// @brief      Gets a pointer to one of the nodes
    //------------------------------------------------------------------------------------------------------- 
    CUDA_HD
    const TreeNode* node_ptr(const size_t i) const { return static_cast<const TreeNode*>(&nodes[i]); }
    
    //------------------------------------------------------------------------------------------------------- 
    /// @brief      Gets a pointer to one of the nodes
    //------------------------------------------------------------------------------------------------------- 
    CUDA_HD
    TreeNode* node_ptr(const size_t i) { return &nodes[i]; }
};

}
}
#endif          // PARAHAPLO_TREE_INTERNAL_H
