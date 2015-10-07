// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo comparators
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_COMPARATORS_HPP
#define PARHAPLO_COMPARATORS_HPP

#include "link_container_cpu.hpp"
#include "node.hpp"

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
    using link_container = const LinkContainer;
    // ------------------------------------------------------------------------------------------------------
private:
    const size_t        _ref_node;          //!< The index that is being used to evaluate against
    link_container&     _links;             //!< The links between the nodes
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the reference node and the links used for comparison
    // ------------------------------------------------------------------------------------------------------
    explicit NodeComparator(const size_t ref_node, link_container& links) 
    : _ref_node(ref_node), _links(links) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Does the comparison between nodes a and b using the reference node
    /// @param[in]  a       The first node for comparison
    /// @param[in]  b       The second node for comparison
    /// @return     The node of higher importance
    // ------------------------------------------------------------------------------------------------------
    bool operator()(const Node& a, const Node& b) const 
    {
        if (a.position() == b.position()) return true;
        
        bool   a_found = true, b_found = true;
        size_t a_value = 0   , b_value = 0; 
        
        // Get the value of a -- casting away atomicity because of read only
        auto lower_idx = std::min(static_cast<const size_t>(a.position()), _ref_node);
        auto upper_idx = std::max(static_cast<const size_t>(a.position()), _ref_node);
        
        if (_links.exists(lower_idx, upper_idx)) {
            auto link = _links.at(lower_idx, upper_idx);      
            a_value   = link.value() - std::min(link.homo_weight()  ,
                                                link.hetro_weight() );
        } else { a_found = false; }
        
    
        // Get the value of b -- casting away atomicity because of read only
        lower_idx = std::min(static_cast<const size_t>(b.position()), _ref_node);
        upper_idx = std::max(static_cast<const size_t>(b.position()), _ref_node);
       
        if (_links.exists(lower_idx, upper_idx)) {
            auto link = _links.at(lower_idx, upper_idx);
            b_value   = link.value() - std::min(link.homo_weight()  ,
                                                link.hetro_weight() );
        } else { b_found = false; }
    
        if (a_found && !b_found) return true;
        if (b_found && !a_found) return false;
        
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

