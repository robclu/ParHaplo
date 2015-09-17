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

// ---------------------------------------------- INPUT FILES -----------------------------------------------

constexpr char* input_simple = "data/test_input_file.txt";

// ----------------------------------------- BLOCK TESTS -----------------------------------------------------

BOOST_AUTO_TEST_SUITE( BlockTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateCpuBlock )
{
    haplo::Block<6          ,                   // 6 rows
                 7          ,                   // 7 columns
                 4          ,                   // 4 cores
                 device::CPU> cpu_block(input_simple);
    
    BOOST_CHECK( cpu_block.num_cores()  == 4 );
    BOOST_CHECK( cpu_block.rows()       == 6 );
    BOOST_CHECK( cpu_block.cols()       == 7 );
}

BOOST_AUTO_TEST_CASE( blockFillsDataWithConstructor )
{  
    // Data filled during construction
    haplo::Block<10, 7, 5, device::CPU> cpu_block(input_simple);
    
    BOOST_CHECK( cpu_block.size() == 70 );
    BOOST_CHECK( cpu_block(0, 0) == 0 );
    BOOST_CHECK( cpu_block(0, 1) == 2 );
    BOOST_CHECK( cpu_block(1, 0) == 1 );
    BOOST_CHECK( cpu_block(1, 1) == 0 );
    BOOST_CHECK( cpu_block(2, 2) == 0 );
    BOOST_CHECK( cpu_block(5, 3) == 0 );
    BOOST_CHECK( cpu_block(8, 5) == 2 ); 
}

BOOST_AUTO_TEST_CASE( blockDeterminesReadInfoCorrectly )
{
    // Read info determined during construction
    haplo::Block<10, 7, 6, device::CPU> cpu_block(input_simple);
    
    BOOST_CHECK( cpu_block[0].start()    == 0 );
    BOOST_CHECK( cpu_block[0].end()      == 0 );
    BOOST_CHECK( cpu_block[0].length()   == 1 );
    BOOST_CHECK( cpu_block[4].start()    == 2 );
    BOOST_CHECK( cpu_block[4].end()      == 3 );
    BOOST_CHECK( cpu_block[4].length()   == 2 ); 
    BOOST_CHECK( cpu_block[8].start()    == 3 );
    BOOST_CHECK( cpu_block[8].end()      == 6 );
    BOOST_CHECK( cpu_block[8].length()   == 4 );
} 

BOOST_AUTO_TEST_CASE( canGetUnsplittableSubBlockInfo )
{
    // Constructor determines subblock info 
    haplo::Block<10, 7, 6, device::CPU> cpu_block(input_simple);

    BOOST_CHECK( cpu_block.subblock_info(0).start() == 0 );    
    BOOST_CHECK( cpu_block.subblock_info(0).end()   == 2 );    
    BOOST_CHECK( cpu_block.subblock_info(1).start() == 2 );    
    BOOST_CHECK( cpu_block.subblock_info(1).end()   == 3 );    
    BOOST_CHECK( cpu_block.subblock_info(2).start() == 3 );    
    BOOST_CHECK( cpu_block.subblock_info(2).end()   == 6 );    
}

BOOST_AUTO_TEST_SUITE_END()

