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
    using iterator          = SearchNode*;
    // ------------------------------------------------------------------------------------------------------
private:
    atomic_type     _free_index;    //!< The index of the next free node
    node_container  _nodes;         //!< The nodes to manage
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    explicit NodeManager() noexcept
    : _free_index{0}, _nodes(0) {}    
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the intiial guess at the number of nodes
    // ------------------------------------------------------------------------------------------------------
    explicit NodeManager(const size_t num_nodes) noexcept
    : _free_index{0}, _nodes(num_nodes) {}

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
    /// @brief      Gets the index of a free node from the container and leaves the following two nodes free 
    ///             for its left and right sub-nodes
    // ------------------------------------------------------------------------------------------------------
    iterator get_new_node();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a reference to a node
    /// @param[in]  index   The index of the element to get
    /// @return     A constant reference to the node
    // ------------------------------------------------------------------------------------------------------
    const_reference node(const size_t index) const;
};


// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

SearchNode* NodeManager<devices::cpu>::get_new_node() 
{
    // Check if we have enough memory 
    if (_free_index + 3 >= _nodes.size()) {}
        // Need to allocate some memory
        
    return &_nodes[_free_index.fetch_and_add(3)];
}

const SearchNode& NodeManager<devices::cpu>::node(const size_t index) const
{
    return _nodes[index];
}


}                   // End namespace haplo
#endif              // PARHAPLO_SEARCH_NODE_HPP