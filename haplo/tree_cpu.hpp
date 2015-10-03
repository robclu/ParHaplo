// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPUHPP

A#include "node_container_cpu.hpp"

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Tree    
/// @brief      Holds nodes which can then be searched to find the optimal haplotypes
/// @tparam     DeviceType  The type of device to use the node on -- so that we can optimize functions for the
///             different implementations and so that each one can be excluded from compilation if necessary
// ----------------------------------------------------------------------------------------------------------
template <>
class Tree<devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
private:
    NodeContainer<devices::cpu>     _nodes;     //!< The nodes in the tree
public:
    Tree();
};

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_HPP

