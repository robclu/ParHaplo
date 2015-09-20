// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class - this does the branching and bounding when solving the
///         haplotype 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPU_HPP

#include "devices.hpp"
#include "node.hpp"

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Tree 
/// @brief      Defines a tree which is searched for the optimal haplotype of an unsplittable block --
///             the implementation method (hardware device which is used) is selected at compile time 
/// @tparam     DeviceType  The type fo device to use -- CPU, GPU, PHI
// ----------------------------------------------------------------------------------------------------------
template <haplo::Device DeviceType>
class Tree;

// Specialize for CPUs
template <> class Tree<haplo::Device::CPU> {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using node_container = std::vector<Node>;
    // ------------------------------------------------------------------------------------------------------
    node_container  _nodes;                 //!< Nodes for this tree
    short           _score;                 //!< The score of the tree -- how good of a solution it is
public:
    Tree(short index) 
    {
        std::cout << "CPU tree\n";
    }
};

}           // End namespace haplo
    
#endif      // PARAHAPLO_TREE_HPP
