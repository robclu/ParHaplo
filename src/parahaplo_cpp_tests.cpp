// ----------------------------------------------------------------------------------------------------------
/// @file   parahaplo_cpp_tests.cpp
/// @brief  Test suites for the c++ implementation of the parahaplo library using Bost.Unit.
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE       ParahaploCppTests
#include <boost/test/unit_test.hpp>

#include "cpp/block.hpp"

#include <iostream>

BOOST_AUTO_TEST_SUITE( BlockTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlockOfAnyTypeAnSize )
{
    phap::Block<float, 3, 4>    block_34;
    phap::Block<int  , 9, 9>    block_99;
    
    BOOST_CHECK( block_34.size() == 12 );
    BOOST_CHECK( block_99.size() == 81 );
}

BOOST_AUTO_TEST_CASE( canAssignAndGetBlockData )
{
    phap::Block<int, 2, 2>      block_22;
    int                         data[4] = {1, 2, 3, 4};
    
    // Assign data - static assert for dimensions mismatch
    block_22.assign_data(data);
    
    // Get a reference to the data
    const std::vector<int>& block_22_data = block_22.get_data();
    
    BOOST_CHECK( block_22_data[0] == 1 );
    BOOST_CHECK( block_22_data[1] == 2 );
    BOOST_CHECK( block_22_data[2] == 3 );
    BOOST_CHECK( block_22_data[3] == 4 );
}

BOOST_AUTO_TEST_SUITE_END()
