// ----------------------------------------------------------------------------------------------------------
/// @file   sorter_tests.cpp
/// @brief  Test suite for parahaplo sorter tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SorterTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/node_container_cpu.hpp"
#include "../haplo/sorter_cpu.hpp"

BOOST_AUTO_TEST_SUITE( SorterSuite )
    
BOOST_AUTO_TEST_CASE( canCorrectlySortNodes )
{
    using link_container = std::vector<haplo::Link>;
    
    // Create a node container with 4 nodes
    haplo::NodeContainer<haplo::devices::cpu> nodes(4);
    
    // Set the weights
    for (size_t i = 0; i < 4; ++i) nodes.weight(i) = i + 2;
   
    // Set the links, they will look like
    // 01 : 1 0 
    // 02 : 2 0
    // 03 : 3 0
    // 12 : 3 1
    // 13 : 4 1
    // 23 : 5 1
    for (size_t i = 0; i < 4; ++i) {
        for(size_t j = i + 1; j < 4; ++j) {
            nodes.link(i, j).homo_weight()  = i + j;
            nodes.link(i, j).hetro_weight() = i;
        }
    }
    
    // Create a comparator   
    // use node 0 as the reference node
    haplo::NodeComparator<link_container> comparator(0, 4, nodes.links());
    
    // Create the sorter
    haplo::Sorter<haplo::devices::cpu> sorter;
    
    sorter(nodes.begin(), nodes.end(), comparator);
    
    for (auto& node : nodes.nodes()) std::cout << node.position() << "\n"; 
}

BOOST_AUTO_TEST_SUITE_END()
