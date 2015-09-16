// ----------------------------------------------------------------------------------------------------------
/// @file   parahaplo_tests.cpp
/// @brief  Test suites for the parahaplo library using Bost.Unit.
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE       ParahaploTests
#include <boost/test/unit_test.hpp>

#include "block.hpp"

#include <iostream>

using device = haplo::Device;

// ----------------------------------------- BLOCK TESTS -----------------------------------------------------

BOOST_AUTO_TEST_SUITE( BlockTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateCpuBlock )
{
    haplo::Block<6          ,                   // 6 rows
                 7          ,                   // 7 columns
                 4          ,                   // 4 cores
                 device::CPU> cpu_block;
    
    BOOST_CHECK( cpu_block.num_cores() == 4 );
}

BOOST_AUTO_TEST_CASE( canFillCpuBlock )
{
    haplo::Block<6, 6, 3, device::CPU> cpu_block;
    //cpu_block.fill(input_file);
    
    BOOST_CHECK( cpu_block.size() == 36 );
}

BOOST_AUTO_TEST_SUITE_END()

