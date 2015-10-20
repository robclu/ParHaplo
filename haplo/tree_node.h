// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node for the gpu search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_H
#define PARHAPLO_TREE_NODE_H

#include <stdint.h>

namespace haplo {
    
struct ALIGN(8) TreeNode {
    
    size_t          haplo_idx   ;
    size_t          root_idx    ; 
    size_t          node_idx    ;
    size_t          alignments  ;
    size_t          pruned      ;
    unsigned int    ubound      ;
    unsigned int    lbound      ;
    unsigned int    min_ubound  ;
    uint16_t        value       ;   // 16 bits for memory alignment 
    uint16_t        prune       ;   // 16 bits for memory alignment
    size_t*         indices     ;
    size_t*         read_ids    ;
    uint8_t*        read_values ;
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constuctor 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    TreeNode() noexcept 
    : haplo_idx{0} , root_idx{0}, node_idx{0}, alignments{0}, pruned{0}     , ubound{0}, lbound{0}, 
      min_ubound{0}, value{0}   , prune{0}   , indices(NULL), read_ids(NULL), read_values(NULL) {} 

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD 
    ~TreeNode() noexcept
    {
        if (alignments > 0) {
            free(indices); free(read_ids); free(read_values);
        }
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Assignment operator
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    void operator=(const TreeNode& other)
    {
        haplo_idx  = other.haplo_idx    ;
        root_idx   = other.root_idx     ;
        alignments = other.alignments   ;
        pruned     = other.pruned       ;
        ubound     = other.ubound       ;
        lbound     = other.lbound       ;
        min_ubound = other.min_ubound   ;
        value      = other.value        ;
        prune      = other.prune        ;
        
        // Pointer copy
        if (alignments != other.alignments) {
            free(indices); free(read_ids); free(read_values); 
          
            indices       = new size_t[other.alignments];
            read_ids      = new size_t[other.alignments];
            read_values   = new uint8_t[other.alignments];
        }
        
        for (size_t i = 0; i < other.alignments; ++i) {
            indices[i]     = other.indices[i];
            read_ids[i]    = other.read_ids[i];
            read_values[i] = other.read_values[i];
        }
    }
};

}

#endif 
