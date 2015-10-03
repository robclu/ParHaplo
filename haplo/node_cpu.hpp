// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo node cpu class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_CPU_HPP
#define PARHAPLO_NODE_CPU_HPP

#include "devices.hpp"

#include <vector>

namespace haplo {
   
// Specialization for the CPU 
template <>
class Node<devices::cpu> {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using link_container    = std::vector<size_t>;      //!< The strength of the links between nodes
private:
    size_t          _weight;        //!< The weight of the node (it's priority)
    size_t          _index;         //!< The haplotype index the node represents
    link_container  _links;         //!< The links between this node and the others
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for when the weight is not known 
    // ------------------------------------------------------------------------------------------------------
    Node() : _weight(1), _index(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor when the weight is known
    /// @param[in]  weight      The weight of the node
    /// @param[in]  index       The index of the node (the position in the haplotype which it represents)
    /// @param[in]  num_nodes   The total number of nodes
    // ------------------------------------------------------------------------------------------------------
    explicit Node(const size_t weight, const size_t index, const size_t num_nodes);
private:
    
};

// ---------------------------------- IMPLEMENTATIONS -------------------------------------------------------

Node<devices::cpu>::Node(const size_t weight, const size_t index, const size_t num_nodes)
: _weight(weight), _index(index), _links(num_nodes - index) 
{
}

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_CPU_HPP

