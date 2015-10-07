// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPU_HPP

#include "devices.hpp"
#include "node_selector_cpu.hpp"
#include "search_node.hpp"
#include "tree.hpp"

#include <iostream>

namespace haplo {
namespace links {
    
static constexpr uint8_t homo   = 0x00;
static constexpr uint8_t hetro  = 0x01;

}               // End namespace links

namespace cores {

static constexpr uint8_t branches   = 0x00;
static constexpr uint8_t opps       = 0x01;
static constexpr uint8_t mixed      = 0x03;

}               // End namespace cores

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
    using link_container    = LinkContainer<devices::cpu>;              // Container for the links
    using atomic_type       = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    atomic_type         _start_node;                //!< The node at which to start the search
    atomic_type         _start_node_worst_case;     //!< The worst case value of the start node
    node_container      _nodes;                     //!< The nodes in the tree
    link_container      _links;                     //!< Links between the nodes of the tree
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    Tree() noexcept : _start_node(0), _start_node_worst_case(0), _nodes(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for a tree
    /// @param[in]  nodes   The number of nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    Tree(const size_t nodes) noexcept : _start_node(0), _start_node_worst_case(0), _nodes(nodes) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Desctructor
    // ------------------------------------------------------------------------------------------------------
    ~Tree() noexcept {}
  
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The mazimum worst case value for the tree
    /// @return     A reference to the maximim worst case value 
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& max_worst_case() { return _start_node_worst_case; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      The index of the start node
    /// @return     A reference to the start node index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& start_node() { return _start_node; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the tree (the number of nodes
    /// @return     The size of of the tree
    // ------------------------------------------------------------------------------------------------------
    inline size_t size() const { return _nodes.num_nodes(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resizes the tree to a certain number of nodes
    /// @param[in]  num_nodes   The number of nodes to create for the tree
    // ------------------------------------------------------------------------------------------------------
    inline void resize(const size_t num_nodes) 
    {
        if (_nodes.num_nodes() != num_nodes) _nodes.resize(num_nodes);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the nodes of the tree
    /// @return     The nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    inline const node_container& nodes() const { return _nodes; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Creates a link for the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)    
    // ------------------------------------------------------------------------------------------------------
    inline void create_link(const size_t node_idx_lower, const size_t node_idx_upper)
    {
        _links.insert(node_idx_lower, node_idx_upper);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the links of the tree
    /// @return     The links for  the tree
    // ------------------------------------------------------------------------------------------------------
    inline const link_container& links() const { return _links; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    /// @tparam     LinkType          The type of the link to get
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t LinkType>
    inline atomic_type& link(const size_t node_idx_lower, const size_t node_idx_upper);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------    
    inline Link& link(const size_t node_idx_lower, const size_t node_idx_upper) 
    {
        return _links.at(node_idx_lower, node_idx_upper);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the max of a link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------
    inline size_t link_max(const size_t node_idx_lower, const size_t node_idx_upper)
    {
        size_t max = 0;
        if (_links.exists(node_idx_lower, node_idx_upper)) 
            max = _links.at(node_idx_lower, node_idx_upper).value();
        else 
            max = 0;
        return max;
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of a node 
    /// @param[in]  idx     The index of the node
    /// @return     The weight of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_weight(const size_t idx) { return _nodes.weight(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the worst case value of a node
    /// @param[in]  idx     The index of the node
    /// @return     The worst case value of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_worst_case(const size_t idx) { return _nodes.worst_case_value(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the haplotype position of a node -- the position in the haplotype a node represents
    /// @param[in]  node_idx    The index of the node to get the haplotype position of
    /// @return     The position the node represents in the haplotype
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_haplo_pos(const size_t node_idx) { return _nodes.haplo_pos(node_idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of a node 
    /// @param[in]  idx     The index of the node
    /// @return     The weight of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& node_weight(const size_t idx) const { return _nodes.weight(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the worst case value of a node
    /// @param[in]  idx     The index of the node
    /// @return     The worst case value of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& node_worst_case(const size_t idx) const { return _nodes.worst_case_value(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the haplotype position of a node -- the position in the haplotype a node represents
    /// @param[in]  node_idx    The index of the node to get the haplotype position of
    /// @return     The position the node represents in the haplotype
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& node_haplo_pos(const size_t node_idx) const { return _nodes.haplo_pos(node_idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Searches the tree for the best solution 
    /// @tparam     CoreMapping     How the cores are mapped for the explaoration -- parallellism is used for
    ///             exploring separate branches, for operation computattions, or a mixture
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t CoreMapping>
    void explore();
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Branches to explore two nodes  
};

// -------------------------------------- IMPLEMENTATIONS ---------------------------------------------------

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::homo>(const size_t node_idx_lower, 
                                                                  const size_t node_idx_upper)
{
    return _links.at(node_idx_lower, node_idx_upper).homo_weight();
}

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::hetro>(const size_t node_idx_lower, 
                                                                   const size_t node_idx_upper)
{
    return _links.at(node_idx_lower, node_idx_upper).hetro_weight();
}

// Not specializing this yet -- going to implement branch based parallelism
template <uint8_t CoreMapping>
void Tree<devices::cpu>::explore() 
{
    // DEBUGGING for the moment
    std::cout << " - - - - - - - EXPLORING TREE - - - - - - -\n";
    
    // The node that's the current reference (for determining how to select nodes)
    atomic_type ref_node = _start_node;
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_CPU_HPP

