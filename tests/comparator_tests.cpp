// ----------------------------------------------------------------------------------------------------------
/// @file   comparator_tests.cpp
/// @brief  Test suite for parahaplo comparator tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE ComparatorTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/comparators.hpp"
#include "../haplo/link_container_cpu.hpp"
#include "../haplo/node_container_cpu.hpp"

BOOST_AUTO_TEST_SUITE( ComparatorSuite )
    
BOOST_AUTO_TEST_CASE( canCorrectlyCompareNodes )
{
    using link_container = haplo::LinkContainer<haplo::devices::cpu>;
    link_container links;
    
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
    for (size_t i = 0; i < 4; ++i) 
        for(size_t j = i + 1; j < 4; ++j) 
            links.insert(i, j, haplo::Link(i + j, i));
    
    // Create a comparator -- use node 0 as the reference node
    haplo::NodeComparator<link_container> comparator(0, links);
    
    BOOST_CHECK( comparator(nodes[1], nodes[2]) == false );
    BOOST_CHECK( comparator(nodes[2], nodes[1]) == true  );
    BOOST_CHECK( comparator(nodes[1], nodes[3]) == false );
    BOOST_CHECK( comparator(nodes[3], nodes[1]) == true  );
    BOOST_CHECK( comparator(nodes[2], nodes[3]) == false );
    BOOST_CHECK( comparator(nodes[3], nodes[2]) == true  );
    
}

BOOST_AUTO_TEST_SUITE_END()
