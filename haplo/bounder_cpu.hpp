// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo bound calculator class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_BOUNDER_CPU_HPP
#define PARHAPLO_BOUNDER_CPU_HPP

#include "bounder.hpp"
#include "link_container_cpu.hpp"
#include "node_container_cpu.hpp"
#include "operations.hpp"
#include "search_node.hpp"
#include <tbb/tbb.h>
#include <algorithm>

namespace haplo {

// Update atomic varibale to max
template <typename T1, typename T2>
void atomic_max_update(tbb::atomic<T1>& atomic_var, T2 value)
{
    T1 state;
    do {
        state = atomic_var;         // Capture state
        if (state >= value) break;  // Exit earlt
    } while (atomic_var.compare_and_swap(value, state) != state);
}

// Specialization for cpu upper bounder
template <>
class Bounder<devices::cpu> {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using link_container    = LinkContainer<devices::cpu>;
    using node_container    = NodeContainer<devices::cpu>;  
    using search_container  = tbb::concurrent_vector<SearchNode>;
    using atomic_type       = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    const node_container&   _nodes;             //!< The nodes used for selection
    const link_container&   _links;             //!< The links used for sorting
    const search_container& _searched_nodes;    //!< Values of nodes which have already been searched
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the parmeters of the bounder
    /// @param[in]  nodes           The nodes used to for calculating the bound
    /// @param[in]  links           The links between the nodes
    /// @param[in]  searched_nodes  The nodes which have already been searched
    // ------------------------------------------------------------------------------------------------------
    Bounder(const node_container& nodes, const link_container& links, const search_container& searched_nodes)
    : _nodes(nodes), _links(links), _searched_nodes(searched_nodes) {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Function which finds the bounds 
    /// @param[in]  node            The node to determine the bound for
    /// @param[in]  node_idx        The index of the node to find the bound for
    /// @param[in]  last_node_idx   The index of the last node to compare with
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t Cores>
    Bounds calculate(const SearchNode& node, const size_t node_idx, const size_t last_node_idx) const;
private: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Finds the value of the comparison node 
    /// @param[in]  thread_idx  The index of the thread doing the searching
    /// @param[in]  start_index The index in the node array to start from
    // ------------------------------------------------------------------------------------------------------
    uint8_t comparison_node_value(const size_t thread_idx, const size_t start_idx) const;
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

template <uint8_t Cores>
Bounds Bounder<devices::cpu>::calculate(const SearchNode& node    ,
                                        const size_t node_idx     ,   
                                        const size_t last_node_idx) const
{
    atomic_type ubound{0};
    atomic_type lbound{0};
    
    const size_t threads = Cores > last_node_idx ? last_node_idx : Cores;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads),
        [&](const tbb::blocked_range<size_t>& thread_ids)
        {
            for (size_t thread_id = thread_ids.begin(); thread_id != thread_ids.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, last_node_idx, threads);
               
                // Iterate over the set nodes to determine the lower and upper bounds  
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t offset      = (it * threads + thread_id);  // Offset for comparison node 
                    size_t other_idx   = _nodes[offset].position();                   // Haplo index of other node
                    size_t lower_bound = 0, upper_bound = 0;

                    // Check if there is a link between the two nodes
                    if (_links.exists(node_idx, other_idx)) {
                       if (comparison_node_value(offset, node.root()) == 0) {   // Other node has a value of 0
                            lower_bound = node.type() == 0 
                                                      ? _links.at(node_idx, other_idx).hetro_weight()
                                                      : _links.at(node_idx, other_idx).homo_weight();
                            upper_bound = node.type() == 0
                                                      ? _links.at(node_idx, other_idx).homo_weight()
                                                      : _links.at(node_idx, other_idx).hetro_weight();
                        } else {                                            // Other node has a value of 1
                            lower_bound = node.type() == 0 
                                                      ? _links.at(node_idx, other_idx).homo_weight()
                                                      : _links.at(node_idx, other_idx).hetro_weight();                  
                            upper_bound = node.type() == 0
                                                      ? _links.at(node_idx, other_idx).hetro_weight()
                                                      : _links.at(node_idx, other_idx).homo_weight(); 
                        }
              
                        // Update the variable 
                        atomic_max_update(lbound, lower_bound);
                        if (lbound == lower_bound) atomic_max_update(ubound, upper_bound);
                    }
                }
            }
        }
    );
    return Bounds(lbound, ubound);
}

uint8_t Bounder<devices::cpu>::comparison_node_value(const size_t thread_idx, const size_t start_idx) const
{
    size_t comparison_idx = start_idx, stop_condition = 0;
    while (stop_condition++ < thread_idx) comparison_idx = _searched_nodes[comparison_idx].root();
    return _searched_nodes[comparison_idx].type();
}

}               // End namespace haplo
#endif          // PARHAPLO_UPPER_BOUNDER_CPU_HPP
