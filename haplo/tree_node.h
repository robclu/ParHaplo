// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node for the gpu search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_H
#define PARHAPLO_TREE_NODE_H

#include <stdint.h>

namespace haplo {
    
struct TreeNode {
    size_t      haplo_idx   ;
    size_t      root_idx    ; 
    size_t      node_idx    ;
    size_t      alignments  ;
    uint8_t     value       ;
    size_t*     read_ids    ;
    uint8_t*    read_values ;
};

}

#endif 
