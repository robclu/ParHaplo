// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node for the gpu search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_H
#define PARHAPLO_TREE_NODE_H

#include <stdint.h>

namespace haplo {
    
// 32 bit stuct, so aligns optimally -- 2 loads for the struct
struct ALIGN(16) TreeNode {
    
    size_t          haplo_idx       ;       // Index of the haplotype the node represents
    size_t          root_idx        ;       // Inddex on the array of nodes of the root node
    size_t          align_idx       ;       // Start index in the array of alignments for this node
    uint32_t        ubound          ;       // The upper bound for the node
    uint32_t        lbound : 31     ;       // The lower bound for the node
    uint32_t        value  : 1      ;       // Value of the node (0,1), 16 bits for memory alignment 
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constuctor 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    TreeNode() 
    : haplo_idx{0} , root_idx{0}, align_idx{0}, ubound{0}, lbound{0}, value{0} {} 

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
        ubound     = other.ubound       ;
        lbound     = other.lbound       ;
        value      = other.value        ;
    }
};

}

#endif 
