// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo search node manager class -- creates an array of search nodes and
///         manages the  nodes
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_MANAGER_HPP
#define PARHAPLO_NODE_MANAGER_HPP

#include "devices.hpp"
#include "search_node.hpp"
#include <tbb/concurrent_vector.h>

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @brief      Manages the search nodes in a tree -- it is efficient search wits, but not efficient cache
///             wise -- room for improvement
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class NodeManager;

// Specialze for the cpu
template <>
class NodeManager<devices::cpu> {
public: 
    // ------------------------------------- ALIAS'S --------------------------------------------------------
    using node_container    = tbb::concurrent_vector<SearchNode>;
    using atomic_type       = tbb::atomic<size_t>;
    using const_reference   = typename node_container::const_reference;
    using reference         = typename node_container::reference;
    // ------------------------------------------------------------------------------------------------------
private:
    atomic_type     _next_index;            //!< The index of the next free node
    node_container  _nodes;                 //!< The nodes to manage
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor -- min 4 nodes
    // ------------------------------------------------------------------------------------------------------
    explicit NodeManager() noexcept
    : _next_index{3}, _nodes(4) {}    
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the intiial guess at the number of nodes
    // ------------------------------------------------------------------------------------------------------
    explicit NodeManager(const size_t num_nodes) noexcept
    : _next_index{3}, _nodes(num_nodes > 4 ? num_nodes : 4) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the node manager -- number of nodes
    /// @return     The number of nodes being managed
    // ------------------------------------------------------------------------------------------------------
    size_t num_nodes() const { return _nodes.size(); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resizes the node container 
    /// @param[in]  num_nodes   The number of nodes to resize to
    // ------------------------------------------------------------------------------------------------------
    void resize(const size_t num_nodes) 
    {
        size_t elements = _nodes.size();
        _nodes.reserve(num_nodes);
        while (elements++ < num_nodes) _nodes.push_back(SearchNode());
    }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a constant reference to the nodes that are being managed
    // ------------------------------------------------------------------------------------------------------
    const node_container& nodes() const { return _nodes; }
     
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the next node, adn makes space for another after it -- the left and
    ///             subnodes of the tree
    /// @return     The index of the next available subnode
    // ------------------------------------------------------------------------------------------------------
    inline size_t get_next_node() { return _next_index.fetch_and_add(2); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a reference to a node
    /// @param[in]  index   The index of the element to get
    /// @return     A constant reference to the node
    // ------------------------------------------------------------------------------------------------------
    inline const_reference node(const size_t index) const { return _nodes[index]; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a reference to a node
    /// @param[in]  index   The index of the element to get
    /// @return     A constant reference to the node
    // ------------------------------------------------------------------------------------------------------
    inline reference node(const size_t index) { return _nodes[index]; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Checks to see if there is enough space for another level of nodes -- alllocates space if
    ///             there is not 
    /// @param[in]  nodes_in_level      The number of nodes in the level
    /// @param[in]  additional_levels   The expected number of additional levels (a guess)
    // ------------------------------------------------------------------------------------------------------
    void add_node_level(const size_t nodes_in_level, const size_t additional_levels = 1);
};


// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

void NodeManager<devices::cpu>::add_node_level(const size_t nodes_in_level, const size_t additional_levels)
{
    // If we don't have enough nodes in the continer
    if (_next_index + (nodes_in_level + 2 * nodes_in_level) >= _nodes.size()) {
        auto additional_level_nodes = [=]() {
            size_t additional_nodes = 2 * nodes_in_level;
            for (size_t i = 1; i < additional_levels; ++i) 
                additional_nodes += 2 * additional_nodes;
            return nodes_in_level + additional_nodes;
        };
 
        size_t total_additional_nodes = 
            _next_index + nodes_in_level - _nodes.size() + additional_level_nodes();
        
        // Create some space and fille with empty nodes 
        _nodes.reserve(_nodes.size() + total_additional_nodes);
        for (size_t i = _next_index; i < total_additional_nodes; ++i) _nodes.push_back(SearchNode());
    }
}



}                   // End namespace haplo
#endif              // PARHAPLO_SEARCH_NODE_HPP