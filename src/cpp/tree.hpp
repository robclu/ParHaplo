// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class - this does the branching and bounding 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_HPP
#define PARHAPLO_TREE_HPP

#include "devices.hpp"

#include <iostream>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Tree 
/// @brief      Defines a tree which is searched for the optimal haplotype of an unsplittable block --
///             the implementation method (hardware device which is used) is selected at compile time 
/// @tparam     DeviceType  The type fo device to use -- CPU, GPU, PHI
// ----------------------------------------------------------------------------------------------------------
template <haplo::Device DeviceType>
class Tree {
};

// Specialize for CPUs
template <> class Tree<haplo::Device::CPU> {
public:
    Tree() 
    {
        std::cout << "CPU tree\n";
    }
};

}           // End namespace haplo
    
#endif      // PARAHAPLO_TREE_HPP
