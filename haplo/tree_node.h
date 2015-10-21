// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node for the gpu search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_H
#define PARHAPLO_TREE_NODE_H

#include <stdint.h>

namespace haplo {
    
// 48 bit stuct, so aligns optimally -- 3 loads for the struct
struct ALIGN(16) TreeNode {
    
    size_t          haplo_idx       ;       // Index of the haplotype the node represents
    size_t          root_idx        ;       // Inddex on the array of nodes of the root node
    size_t          pruned          ;       // Number of pruned nodes based on this node's lower bound
    size_t          align_idx       ;       // Start index in the array of alignments for this node
    unsigned int    ubound          ;       // The upper bound for the node
    unsigned int    lbound          ;       // The lower bound for the node
    unsigned int    min_ubound      ;       // The minimum lower bound seen by the node
    uint16_t        value : 1       ;       // Value of the node (0,1), 16 bits for memory alignment 
    uint16_t        prune : 1       ;       // If the node is pruned, 16 bits for memory alignment
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constuctor 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    TreeNode() noexcept 
    : haplo_idx{0} , root_idx{0}, align_idx{0} , pruned{0}, 
      ubound{0}    , lbound{0}  , min_ubound{0}, value{0} , prune{0}   {} 

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD 
    ~TreeNode() noexcept {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Assignment operator
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    void operator=(const TreeNode& other)
    {
        haplo_idx  = other.haplo_idx    ;
        root_idx   = other.root_idx     ;
        pruned     = other.pruned       ;
        ubound     = other.ubound       ;
        lbound     = other.lbound       ;
        min_ubound = other.min_ubound   ;
        value      = other.value        ;
        prune      = other.prune        ;
    }
};

}

#endif 
