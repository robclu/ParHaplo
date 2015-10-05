// ----------------------------------------------------------------------------------------------------------
/// @file   tree_tests.cpp
/// @brief  Test suite for parahaplo tree tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE TreeTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/tree_cpu.hpp"

BOOST_AUTO_TEST_SUITE( TreeSuite )
    
BOOST_AUTO_TEST_CASE( canCreateATreeAndGetAndSetLinks )
{
    // Create a tree with 12 nodes
    haplo::Tree<haplo::devices::cpu> tree(12);
 
    tree.link<haplo::links::homo>(0, 1)  = 5;
    tree.link<haplo::links::hetro>(0, 1) = 4;
    
    tree.link<haplo::links::hetro>(0, 3) += 4;
    
    BOOST_CHECK( tree.link<haplo::links::homo>(0, 1)   == 5 );
    BOOST_CHECK( tree.link<haplo::links::homo>(0, 2)   == 0 );
    BOOST_CHECK( tree.link<haplo::links::hetro>(0, 1)  == 4 );
    BOOST_CHECK( tree.link<haplo::links::hetro>(0, 2)  == 0 );
    BOOST_CHECK( tree.link<haplo::links::hetro>(0, 3)  == 4 );
    BOOST_CHECK( tree.link_max(0, 1)                   == 5 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodeHaploTypePosition )
{
    // Create a tree with 12 nodes
    haplo::Tree<haplo::devices::cpu> tree(12);    
    
    BOOST_CHECK( tree.node_haplo_pos(4) == 4 );
    
    tree.node_haplo_pos(4) = 12;
    
    BOOST_CHECK( tree.node_haplo_pos(4) == 12 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodeWeight )
{
    // Create a tree with 12 nodes
    haplo::Tree<haplo::devices::cpu> tree(12);    
    
    BOOST_CHECK( tree.node_weight(4) == 1 );
    
    tree.node_weight(4) = 12;
    
    BOOST_CHECK( tree.node_weight(4) == 12 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodeWorstCaseValue )
{
    // Create a tree with 12 nodes
    haplo::Tree<haplo::devices::cpu> tree(12);    
    
    tree.node_worst_case(2) = 12;
    
    BOOST_CHECK( tree.node_worst_case(2) == 12 );    
    BOOST_CHECK( tree.node_worst_case(4) == 0  );    
}

BOOST_AUTO_TEST_SUITE_END()
