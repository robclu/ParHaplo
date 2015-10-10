// ----------------------------------------------------------------------------------------------------------
/// @file   search_node_tests.cpp
/// @brief  Test suite for parahaplo search node tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SearchNodeTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/search_node.hpp"

BOOST_AUTO_TEST_SUITE( SearchNodeSuite )

BOOST_AUTO_TEST_CASE( canGetIndex ) 
{
    haplo::SearchNode node(12, 0, 1, 1111, 8394);
    
    BOOST_CHECK( node.index() == 12 );
}

BOOST_AUTO_TEST_CASE( canGetValue ) 
{
    haplo::SearchNode node(12, 0, 1, 1111, 8394);
    
    BOOST_CHECK( node.value() == 0 );
}

BOOST_AUTO_TEST_CASE( canHandleValueOutOfRange1 ) 
{
    // Value = 3 is too big, must be binary
    haplo::SearchNode node(12, 3, 0, 1111, 8394);
    
    // Should set value == 1 since it ands with 0x01
    BOOST_CHECK( node.value() == 1 );
}

BOOST_AUTO_TEST_CASE( canHandleValueOutOfRange2 ) 
{
    // Value = 4 is too big, must be binary
    haplo::SearchNode node(12, 4, 1, 1111, 8394);
    
    // Should set value == 0 since it ands with 0x01
    BOOST_CHECK( node.value() == 0 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetUpperBound ) 
{
    haplo::SearchNode node(12, 0, 1, 64998, 8394);
    
    BOOST_CHECK( node.upper_bound() == 64998 );
    
    node.upper_bound() = 77;
    
    BOOST_CHECK( node.upper_bound() == 77 );   
}

BOOST_AUTO_TEST_CASE( canGetAndSetLowerBound ) 
{
    haplo::SearchNode node(12, 0, 0, 1111, 65012);
    
    BOOST_CHECK( node.lower_bound() == 65012 );
    
    node.lower_bound() = 7;
    
    BOOST_CHECK( node.lower_bound() == 7 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetLeftNode ) 
{
    haplo::SearchNode node(12, 0, 0, 1111, 65012);
    
    node.left() = 14;
    
    BOOST_CHECK( node.left() == 14 );
    
}

BOOST_AUTO_TEST_CASE( canGetAndSetRightNode ) 
{
    haplo::SearchNode node(12, 0, 1, 1111, 65012);
 
    node.right() = 999;
    
    BOOST_CHECK( node.right() = 999 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetType ) 
{
    haplo::SearchNode node(12, 0, 1, 1111, 65012);
 
    node.set_type(haplo::types::left);
    
    BOOST_CHECK( node.type() == haplo::types::left );
    
    node.set_type(haplo::types::right);
    
    BOOST_CHECK( node.type() == haplo::types::right );
}

BOOST_AUTO_TEST_CASE( canGetAndSetRoot ) 
{
    haplo::SearchNode node;
 
    node.root() = 4;
    
    BOOST_CHECK( node.root() == 4 );
}

BOOST_AUTO_TEST_SUITE_END()
