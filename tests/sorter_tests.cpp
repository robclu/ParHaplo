// ----------------------------------------------------------------------------------------------------------
/// @file   sorter_tests.cpp
/// @brief  Test suite for parahaplo sorter tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SorterTests
#endif
#include <boost/test/unit_test.hpp>
#include <chrono>
#include <iostream>

#include "../haplo/node_container_cpu.hpp"
#include "../haplo/sorter_cpu.hpp"

using namespace std::chrono;

BOOST_AUTO_TEST_SUITE( SorterSuite )
    
BOOST_AUTO_TEST_CASE( canCorrectlySortNodes )
{
    const size_t num_nodes = 350;
    
    using link_container = std::vector<haplo::Link>;
    
    // Create a node container with 4 nodes
    haplo::NodeContainer<haplo::devices::cpu> nodes(num_nodes);
    
    // Set the weights
    for (size_t i = 0; i < num_nodes; ++i) nodes.weight(i) = i + 2;
   
    // Set the links
    for (size_t i = 0; i < num_nodes; ++i) {
        for(size_t j = i + 1; j < num_nodes; ++j) {
            nodes.link(i, j).homo_weight()  = i + j;
            nodes.link(i, j).hetro_weight() = i;
        } 
    }
    
    // Create a comparator   
    // use node 0 as the reference node
    haplo::NodeComparator<link_container> comparator(0, num_nodes, nodes.links());
    
    // Create the sorter
    haplo::Sorter<haplo::devices::cpu> sorter;

    // Time the sort    
    high_resolution_clock::time_point start = high_resolution_clock::now(); 
    
    // Sorted nodes should be in the order of 
    //  0 3 2 1
    sorter(nodes.begin(), nodes.end(), comparator);
    
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> sort_time = duration_cast<duration<double>>(end - start);
   
    size_t node_idx = num_nodes - 1; 
    for (auto& node : nodes.nodes()) {
        if (node_idx != num_nodes - 1) {   // Not comparing the first node
            BOOST_CHECK( node.position() == node_idx );
            --node_idx;
        }
    }
    
    std::cout << "Sort time : " << sort_time.count() << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
