// ----------------------------------------------------------------------------------------------------------
/// @file   node_selector_tests.cpp
/// @brief  Test suite for parahaplo node selector tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE NodeSelectorTests
#endif
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <iostream>

#include "../haplo/node_selector_cpu.hpp"

using namespace std::chrono;

BOOST_AUTO_TEST_SUITE( NodeSelectorSuite )
    
BOOST_AUTO_TEST_CASE( canCorrectlySelectNodes )
{
    const size_t num_nodes = 4;
    
    using link_container = haplo::LinkContainer<haplo::devices::cpu>;
    link_container links;
    
    // Create a node container with 4 nodes
    haplo::NodeContainer<haplo::devices::cpu> nodes(num_nodes);
    
    // Set the weights
    for (size_t i = 0; i < num_nodes; ++i) nodes.weight(i) = i + 2;
   
    // Set the links
    // 01 : 1 0
    // 02 : 2 0
    // 03 : 3 0
    // 12 : 3 1
    // 13 : 4 1
    // 23 : 5 2
    for (size_t i = 0; i < num_nodes; ++i) 
        for(size_t j = i + 1; j < num_nodes; ++j) 
            links.insert(i, j, haplo::Link(i + j, i));

    // Create the selector -- start at node 1
    haplo::NodeSelector<haplo::devices::cpu> selector(nodes, links, 1);

    std::cout << "Starting selection\n";
    
    // Time the selection
    high_resolution_clock::time_point start = high_resolution_clock::now(); 
    
    // Sorted nodes should be in the order of
    for (int i = 1; i <= num_nodes; ++i)
        BOOST_CHECK( selector.select_node() == i );
    
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> sort_time = duration_cast<duration<double>>(end - start);
   
    std::cout << "Selection time : " << sort_time.count() << "\n";
}

BOOST_AUTO_TEST_CASE( canCorrectlySelectNodes1 )
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
    
    // Make node 3 NIH
    nodes[3].type() = 0;
    
    // This should order the nodes as (choosing 1 to be the start)
    //  1 - 4 - 3 - 0 - 5 - 2 (but 3 is NIH) so
    //  1 - 4 - 0 - 5 - 2 - 3

    // Create the selector -- start at node 1
    haplo::NodeSelector<haplo::devices::cpu> selector(nodes, links, 1);

    std::cout << "Starting selection\n";
    
    // Time the selection
    high_resolution_clock::time_point start = high_resolution_clock::now(); 

    // The first time the selection is called it pointes to the second element 
    BOOST_CHECK( nodes[selector.select_node()].position() == 4);
    BOOST_CHECK( nodes[selector.select_node()].position() == 0);
    BOOST_CHECK( nodes[selector.select_node()].position() == 5);
    BOOST_CHECK( nodes[selector.select_node()].position() == 2);
    BOOST_CHECK( nodes[selector.select_node()].position() == 3);
    
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> sort_time = duration_cast<duration<double>>(end - start);
   
    std::cout << "Selection time : " << sort_time.count() << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
