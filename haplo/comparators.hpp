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
        if (a.type() == 2) return false;                                        // A is a duplicateDuplicate
        if (b.type() == 2) return true;                                         // A is a duplicateDuplicate
        if (a.type() == 1 && b.type() == 0) return false;                       // A is an NIH column
        if (a.type() == 0 && b.type() == 1) return true;                        // B is an NIH column
        if (a.position() == b.position()) return true;                          // Same node
        
        bool   a_found = true, b_found = true;
        size_t a_value = 0, b_value = 0, a_min, b_min;
        
        if (_links.exists(a.position(), _ref_node)) {
            auto link = _links.at(a.position(), _ref_node);      
            a_value   = link.max() - link.min();
            a_min     = link.min(); 
        } else { a_found = false; }
       
        if (_links.exists(b.position(), _ref_node)) {
            auto link = _links.at(b.position(), _ref_node);
            b_value   = link.max() - link.min();
            b_min     = link.min();
        } else { b_found = false; }
    
        if (a_found && !b_found) return true;                           // Link for a but not b
        if (b_found && !a_found) return false;                          // Link for b but not a 
        
        return a_value > b_value ? true 
                                 : b_value > a_value 
                                    ? false 
                                    : a_min < b_min 
                                        ? true
                                        : b_min < a_min
                                            ? false 
                                            : a.weight() >= b.weight() 
                                                ? true 
                                                : false;
    }
};

}               // End namespace haplo
#endif          // PARAHAPLO_COMPARATORS_HPP

