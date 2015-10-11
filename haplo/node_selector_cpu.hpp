// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo node selector class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_SELECTOR_CPU_HPP
#define PARHAPLO_NODE_SELECTOR_CPU_HPP

#include "comparators.hpp"
#include "link_container_cpu.hpp"
#include "node_container_cpu.hpp"
#include "node_selector.hpp"
#include "sorter_cpu.hpp"

#include <iostream>

namespace haplo {

// Specialization for cpu selection operator
template <>
class NodeSelector<devices::cpu> {
public:
    // ---------------------------------------- ALIAS'S -----------------------------------------------------
    using sorter_type       = Sorter<devices::cpu>;
    using link_container    = LinkContainer<devices::cpu>;
    using node_container    = NodeContainer<devices::cpu>;
    using node_comparator   = NodeComparator<link_container>;
    // ------------------------------------------------------------------------------------------------------
private:
    sorter_type             _sorter;        //!< The sorter for the nodes
    node_container&         _nodes;         //!< The nodes used for selection<F12>
    const link_container&   _links;         //!< The links used for sorting
    size_t                  _ref_node;      //!< The node being used as the reference
    size_t                  _next_node;     //!< The index of the node to look at next
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the parameters of the node selector
    /// @param[in]  nodes       The nodes available for selection
    /// @param[in]  links       The links used to determine the best selection of the nodes
    /// @param[in]  start_node  The node to start with (used as the root of the tree)
    // ------------------------------------------------------------------------------------------------------
    NodeSelector(node_container& nodes, const link_container& links, const size_t start_node = 0);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Selects the next node, may need to sort the nodes on occasion, and updates the state of
    ///             node selector
    /// @param      node_index  The index of the node in the node container
    // ------------------------------------------------------------------------------------------------------
    size_t select_node();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the position of the last selected node
    /// @rerturn    The index of the last selected node
    // ------------------------------------------------------------------------------------------------------
    inline size_t last_selected_index() const { return _next_node - 1; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of last node to be searched
    /// @rerturn    The index of the last node to search
    // ------------------------------------------------------------------------------------------------------
    inline size_t last_search_index() const { return _nodes.num_nodes(); }

};

// -------------------------------------------- IMPLEMENTATIONS ---------------------------------------------
 
NodeSelector<devices::cpu>::NodeSelector(node_container&        nodes       , 
                                         const link_container&  links       , 
                                         const size_t           start_node  ) 
: _nodes(nodes), _links(links), _ref_node(0), _next_node(1) 
{
    // Sort the nodes based on the start node
    std::iter_swap(nodes.begin(), nodes.begin() + start_node);              // Make the start node the head
    node_comparator comparator(start_node, links);                          // Create a comparator
    _sorter(nodes.begin() + 1, nodes.end(), comparator);                    // Sort the other nodes
    std::cout << "Sorted!\n";
    for (auto node : _nodes.nodes())
        std::cout << node.position() << "\n";
}

 
size_t NodeSelector<devices::cpu>::select_node()
{
    // If next node is not related to ref node -- resort
    if (!_links.exists(_ref_node, _next_node) && _next_node < _nodes.num_nodes()) {
        ++_ref_node;                                                        // Move the referene node up
        node_comparator comparator(_ref_node, _links);                      // Create a comparator
        _sorter(_nodes.begin() + _next_node, _nodes.end(), comparator);     // Sort the unsearched nodes
    }
    // Next node it the index of the next node to be searched
    return _next_node++;
}

}               // End namespace haplo
#endif          // PARHAPLO_NODE_SELECTOR_CPU_HPP

