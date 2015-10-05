// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPU_HPP

#include "devices.hpp"
#include "node_container_cpu.hpp"
#include "tree.hpp"

namespace haplo {
namespace links {
    
static constexpr uint8_t homo   = 0x00;
static constexpr uint8_t hetro  = 0x01;

}               // End namespace links

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
    using node_container    = NodeContainer<devices::cpu>;              // Container for the nodes
    using info_container    = typename node_container::info_container;  // Vector of nodes
    using link_container    = typename node_container::link_container;  // Vector of links
    // ------------------------------------------------------------------------------------------------------
private:
    node_container     _nodes;     //!< The nodes in the tree
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for a tree
    /// @param[in]  nodes   The number of nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    Tree(const size_t nodes) : _nodes(nodes) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the nodes of the tree
    /// @return     The nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    inline const info_container& nodes() const { return _nodes.nodes(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the links of the tree
    /// @return     The links for  the tree
    // ------------------------------------------------------------------------------------------------------
    inline const link_container& links() const { return _nodes.links(); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    /// @tparam     LinkType            Selector for which of the link types to get the value of
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t LinkType>
    inline tbb::atomic<size_t>& link(const size_t node_idx_lower, const size_t node_idx_upper);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the max of a link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------
    inline size_t link_max(const size_t node_idx_lower, const size_t node_idx_upper)
    {
        return _nodes.link(node_idx_lower, node_idx_upper).value();
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of a node 
    /// @param[in]  idx     The index of the node
    /// @return     The weight of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline size_t& node_weight(const size_t idx) { return _nodes.weight(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the haplotype position of a node -- the position in the haplotype a node represents
    /// @param[in]  node_idx    The index of the node to get the haplotype position of
    /// @return     The position the node represents in the haplotype
    // ------------------------------------------------------------------------------------------------------
    inline size_t& node_haplo_pos(const size_t node_idx) { return _nodes.haplo_pos(node_idx); }
};

// -------------------------------------- IMPLEMENTATIONS ---------------------------------------------------

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::homo>(const size_t node_idx_lower, 
                                                                  const size_t node_idx_upper)
{
    return _nodes.link(node_idx_lower, node_idx_upper).homo_weight();
}

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::hetro>(const size_t node_idx_lower, 
                                                                   const size_t node_idx_upper)
{
    return _nodes.link(node_idx_lower, node_idx_upper).hetro_weight();
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_CPU_HPP

