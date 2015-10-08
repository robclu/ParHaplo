// ----------------------------------------------------------------------------------------------------------
/// @file   bounder_tests.cpp
/// @brief  Test suite for parahaplo bound calculation functor tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BounderTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/node_selector_cpu.hpp"
#include "../haplo/bounder_cpu.hpp"

BOOST_AUTO_TEST_SUITE( BounderSuite )

BOOST_AUTO_TEST_CASE( canCorrectlyDetermineBounds )
{
    const size_t num_nodes = 6;
    
    using link_container = haplo::LinkContainer<haplo::devices::cpu>;
    link_container links;
    
    // Create a node container with 4 nodes
    haplo::NodeContainer<haplo::devices::cpu> nodes(num_nodes);
    
    // Set the weights
    for (size_t i = 0; i < num_nodes; ++i) nodes.weight(i) = i + 2;
  
    // Set some links manually 
    links.insert(0, 1, haplo::Link(2, 0));  // 01 : 2 0
    links.insert(1, 3, haplo::Link(5, 0));  // 13 : 5 0
    links.insert(1, 4, haplo::Link(8, 0));  // 14 : 8 0
    links.insert(0, 2, haplo::Link(4, 0));  // 02 : 3 0
    links.insert(0, 4, haplo::Link(0, 3));  // 04 : 0 3
    links.insert(3, 4, haplo::Link(2, 6));  // 34 : 2 6
    
    // This should order the nodes as (choosing 1 to be the start)
    //  1 - 4 - 3 - 0 - 2 - 5

    // Create the selector -- start at node 1
    haplo::NodeSelector<haplo::devices::cpu> selector(nodes, links, 1);

    // Create the bounder object 
    haplo::Bounder<haplo::devices::cpu> bound_calculator(nodes, links);
    
    auto idx        = selector.select_node();
    auto node       = nodes[idx];
    auto end_idx    = selector.last_selected_index();
    
    auto bounds = bound_calculator.calculate<1>(node.position(), end_idx);

    BOOST_CHECK( bounds.lower == 0 );
    BOOST_CHECK( bounds.upper == 8 );
    
    // Select a new node
    idx     = selector.select_node();
    node    = nodes[idx];
    end_idx = selector.last_selected_index();
    
    bounds = bound_calculator.calculate<4>(node.position(), end_idx);
    
    BOOST_CHECK( bounds.lower == 2 );
    BOOST_CHECK( bounds.upper == 6 );
}

BOOST_AUTO_TEST_SUITE_END()
