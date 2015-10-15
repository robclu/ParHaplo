// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for internal tree sruct to clean the GPU interface
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_INTERNAL_H
#define PARHAPLO_TREE_INTERNAL_H

#include "read_info_gpu.cuh"
#include "snp_info_gpu.h"
#include "tree_node_manager.h"

namespace haplo     {
namespace internal  {
 
struct Tree {    
public:
    //-------------------------------------------------------------------------------------------------------
    using read_info_type                = ReadInfo;
    using snp_info_type                 = SnpInfoGpu;
    using small_type                    = uint8_t;
    using node_container                = NodeManagerGpu;
    using small_container               = thrust::host_vector<small_type>;
    //------------------------------------------------------------------------------------------------------- 
   
    //------------------------------------------------------------------------------------------------------- 
    /// @brief  Constructor to initialize avriables
    //------------------------------------------------------------------------------------------------------- 
    Tree(const size_t num_snps, const size_t num_reads)
    : snps(num_snps), reads(num_reads), node_manager(snps), last_searched_snp(0) {}

    size_t                      snps;
    size_t                      reads;
    size_t                      last_searched_snp;
    size_t*                     search_snps;
    small_type*                 data;
    small_type*                 haplotype;
    small_type*                 alignments;
    read_info_type*             read_info;
    snp_info_type*              snp_info;
    node_container              node_manager; 
};

}
}
#endif          // PARAHAPLO_TREE_INTERNAL_H
