// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo search operations class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SEARCH_OPERATIONS_CPU_HPP
#define PARHAPLO_SEARCH_OPERATIONS_CPU_HPP

#include "node_container.hpp"
#include "sorter_cpu.hpp"

namespace haplo {
            

// Specialization for cpu selection operator
template <>
class NodeSelector<devices::cpu> {
public:
    // ---------------------------------------- ALIAS'S -----------------------------------------------------
    using sorter_type       = Sorter<devices::cpu>;
    using links_reference   = std::vector<Link>&;
    // ------------------------------------------------------------------------------------------------------
private:
    sorter_type     _sorter;        //!< The sorter for the nodes
    links_reference _links;         //!< The links used for sorting
    const size_t    _total_nodes;   //!< The total number of nodes 
    size_t          _last_node;     //!< The last selected node
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the parameters of the node selector
    /// @param[in]  total_nodes     The total number of nodes 
    /// @param[in]  links           The links used to determine the best selection of the nodes
    // ------------------------------------------------------------------------------------------------------
    NodeSelector(const size_t total_nodes, link_container& links)
    : _total_nodes(total_nodes), _links(links), _last_node(0) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Selects the next node, may need to sort the nodes on occasion
    /// @param      start_node      The node to start with -- default to 1
    // ------------------------------------------------------------------------------------------------------
    void select_node(const size_t start_node = 0);

};

template <>
void NodeSelector<cpu::devices>::select_node(const size_t start_node)
{
    
}

}               // End namespace haplo
#endif          // PARHAPLO_SEARCH_OPERATIONS_CPU_HPP

