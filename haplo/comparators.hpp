// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo comparators
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_COMPARATORS_HPP
#define PARHAPLO_COMPARATORS_HPP

#include "node_container.hpp"

#include <algorithm>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @brief      Class for comparing nodes 
/// @tparam     LinkContainer   The type of container used to store the links between nodes
// ----------------------------------------------------------------------------------------------------------
template <typename LinkContainer>
class NodeComparator {
public: 
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using links_reference = const LinkContainer&;
    // ------------------------------------------------------------------------------------------------------
private:
    const size_t        _ref_node;          //!< The index that is being used to evaluate against
    const size_t        _nodes;             //!< The number of nodes which made up the links
    links_reference     _links;             //!< The links between the nodes
public:
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    explicit NodeComparator(const size_t ref_node, const size_t nodes, links_reference links) 
    : _ref_node(ref_node), _nodes(nodes), _links(links) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Does the comparison between nodes a and b using the reference node
    /// @param[in]  a       The first node for comparison
    /// @param[in]  b       The second node for comparison
    /// @return     The node of higher importance
    // ------------------------------------------------------------------------------------------------------
    bool operator()(const Node& a, const Node& b) const 
    {
        // Get the value of a -- casting away atomicity because of read only
        auto lower_idx = std::min(static_cast<const size_t>(a.position()), _ref_node);
        auto upper_idx = std::max(static_cast<const size_t>(a.position()), _ref_node);
        
        auto link_idx  = _links.size() - ((_nodes - lower_idx) * (_nodes - lower_idx - 1) / 2)
                           + upper_idx - lower_idx - 1;
        
        auto a_value = _links[link_idx].value() - std::min(_links[link_idx].homo_weight()  ,
                                                           _links[link_idx].hetro_weight() );
    
        // Get the value of b -- casting away atomicity because of read only
        lower_idx = std::min(static_cast<const size_t>(b.position()), _ref_node);
        upper_idx = std::max(static_cast<const size_t>(b.position()), _ref_node);
        
        link_idx  = _links.size() - ((_nodes - lower_idx) * (_nodes - lower_idx - 1) / 2)
                           + upper_idx - lower_idx - 1;
        
        auto b_value = _links[link_idx].value() - std::min(_links[link_idx].homo_weight()  ,
                                                           _links[link_idx].hetro_weight() );
        
        return a_value > b_value ? true 
                                 : b_value > a_value 
                                    ? false 
                                    : a.weight() >= b.weight() 
                                        ? true 
                                        : false;
    }
    
};

}               // End namespace haplo
#endif          // PARAHAPLO_COMPARATORS_HPP

