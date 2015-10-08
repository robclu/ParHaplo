// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo bound calculator class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_BOUNDER_CPU_HPP
#define PARHAPLO_BOUNDER_CPU_HPP

#include "bounder.hpp"
#include "link_container_cpu.hpp"
#include "node_container_cpu.hpp"
#include "operations.hpp"
#include <tbb/tbb.h>
#include <algorithm>

namespace haplo {

// Specialization for cpu upper bounder
template <>
class Bounder<devices::cpu> {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using link_container    = LinkContainer<devices::cpu>;
    using node_container    = NodeContainer<devices::cpu>;  
    using atomic_type       = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    const node_container&   _nodes;             //!< The nodes used for selection
    const link_container&   _links;             //!< The links used for sorting
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the parmeters of the bounder
    /// @param[in]  nodes       The nodes used to for calculating the bound
    /// @param[in]  links       The links between the nodes
    // ------------------------------------------------------------------------------------------------------
    Bounder(const node_container& nodes, const link_container& links) 
    : _nodes(nodes), _links(links) {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Function which finds the bounds 
    /// @param[in]  node_idx        The index of the node to find the bound for
    /// @param[in]  last_node_idx   The index of the last node to compare with
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t Cores>
    Bounds calculate(const size_t node_idx, const size_t last_node_idx) const;
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

template <uint8_t Cores>
Bounds Bounder<devices::cpu>::calculate(const size_t node_idx, const size_t last_node_idx) const
{
    atomic_type max_value{0};
    atomic_type min_value{0};
    
    const size_t threads = Cores > last_node_idx ? last_node_idx : Cores;
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads),
        [&](const tbb::blocked_range<size_t>& thread_ids)
        {
            for (size_t thread_id = thread_ids.begin(); thread_id != thread_ids.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, last_node_idx, threads);
               
                // Iterate over the set nodes to determine the lower and upper bounds  
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t offset      = it * threads + thread_id;              // Offset for comparison node 
                    size_t other_idx   = _nodes[offset].position();             // Haplo index of other node
                    size_t small_idx   = std::min(node_idx, other_idx);         // Smaller of the inddices
                    size_t big_idx     = std::max(node_idx, other_idx);         // Bigger of the indices
                    size_t link_max    = 0, link_min = 0;                       // Max and min value
                    
                    // Check if there is a link between the two nodes
                    if (_links.exists(small_idx, big_idx)) {
                        link_max = _links.at(small_idx, big_idx).value();
                        link_min = _links.at(small_idx, big_idx).min();        
                        
                        // Make sure that if the max is the max, that it is correctly updated
                        size_t state_max;
                        do {
                            state_max = max_value;
                            if (state_max >= link_max) break;
                        } while (max_value.compare_and_swap(link_max, state_max) != state_max);
                        // Set the min value if the max succeeded
                        if (max_value == link_max) min_value.fetch_and_store(link_min);
                    }
                }
            }
        }
    );
    return Bounds(min_value, max_value);
}

}               // End namespace haplo
#endif          // PARHAPLO_UPPER_BOUNDER_CPU_HPP
