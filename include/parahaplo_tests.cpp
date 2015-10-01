// ----------------------------------------------------------------------------------------------------------
/// @file   parahaplo_tests.cpp
/// @brief  Test suites for the parahaplo library using Bost.Unit.
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE       ParahaploTests
#include <boost/test/unit_test.hpp>

#include "block.hpp"
#include "node.hpp"
#include "usplit_block.hpp"
#include "variables.hpp"

#include <iostream>

using device = haplo::Device;

// ---------------------------------------------- INPUT FILES -----------------------------------------------

constexpr char* input_simple    = "data/test_input_file.txt";
constexpr char* input_dup_rows  = "data/test_input_file_duplicate_rows.txt";

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

// ----------------------------------- UNSPLITTABLE BLOCK TESTS ---------------------------------------------

BOOST_AUTO_TEST_SUITE( UnsplittableBlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateUnplittabelBlockAndAccessElements ) 
{
    using block_type = haplo::Block<14, 7, 8, device::CPU>;
    block_type cpu_block(input_dup_rows); 
    
    // By default the unplittable block inherits the device type and the number of 
    // cores to use from the block it uses (as is created below), but we could specify 
    // the device type and number of cores explicitly as template parameters, such as
    // 
    // haplo::UnsplittableBlock<block_type, Device::GPU, 1000>
    // 
    // To use a GPU with 1000 cores
    haplo::UnsplittableBlock<block_type> usplit_block(cpu_block, 0);  
    
    BOOST_CHECK( usplit_block.at(0, 0) == 1 );
    BOOST_CHECK( usplit_block.at(0, 1) == 0 );
    BOOST_CHECK( usplit_block.at(0, 2) == 2 );
    BOOST_CHECK( usplit_block.at(1, 0) == 1 );
    BOOST_CHECK( usplit_block.at(1, 1) == 0 );
    BOOST_CHECK( usplit_block.at(1, 2) == 2 );
    BOOST_CHECK( usplit_block.at(2, 0) == 0 );
    BOOST_CHECK( usplit_block.at(3, 0) == 0 );
    BOOST_CHECK( usplit_block.at(4, 0) == 0 );
    BOOST_CHECK( usplit_block.at(5, 1) == 0 );
}

BOOST_AUTO_TEST_CASE( canCreateAnUnsplittableBlockAndDetermineRowMultiplicities )
{
    using block_type = haplo::Block<14, 7, 8, device::CPU>;
    block_type cpu_block(input_dup_rows); 
    
    haplo::UnsplittableBlock<block_type> usplit_block(cpu_block, 0); 

    BOOST_CHECK( usplit_block.row_multiplicity(0) == 2 );
    BOOST_CHECK( usplit_block.row_multiplicity(1) == 3 );
    BOOST_CHECK( usplit_block.row_multiplicity(2) == 1 );
}

// TODO: Add column tests

BOOST_AUTO_TEST_CASE( canSolveHaplotypeForUnsplittableBlock )
{
    using block_type = haplo::Block<14, 7, 8, device::CPU>;
    block_type cpu_block(input_dup_rows);
    
    haplo::UnsplittableBlock<block_type> ublock(cpu_block, 0);
    
    ublock.solve_haplotype();
}

BOOST_AUTO_TEST_SUITE_END()
    
// ------------------------------------------------- NODE TESTS ---------------------------------------------

BOOST_AUTO_TEST_SUITE( NodeSuite )

BOOST_AUTO_TEST_CASE( canCreateANode )
{
    haplo::Node node_correct1(  0, 0, 0, 0 );
    haplo::Node node_correct2(  0, 0, 1, 1 );
    haplo::Node node_correct3(  0, 0, 0, 1 );
    haplo::Node node_incorrect( 0, 0, 2, 2 );
    
    BOOST_CHECK( node_correct1.x_value() == 0 );
    BOOST_CHECK( node_correct1.y_value() == 0 );
    BOOST_CHECK( node_correct2.x_value() == 1 );
    BOOST_CHECK( node_correct2.y_value() == 1 );
    BOOST_CHECK( node_incorrect.x_value() == 0 );
    BOOST_CHECK( node_incorrect.y_value() == 0 );
}

BOOST_AUTO_TEST_SUITE_END()
    
// -------------------------------------------- VARIABLES TESTS ---------------------------------------------

BOOST_AUTO_TEST_SUITE( VariablesSuite )

BOOST_AUTO_TEST_CASE( canCreateABinaryArrayAndSetElementANode )
{
    haplo::BinaryArray bit_array;
    
    bit_array.set(4, 1); 
    bit_array.set(3, 1); 
    bit_array.set(0, 1); 
    
    BOOST_CHECK( bit_array[0] == 1 );
    BOOST_CHECK( bit_array[1] == 0 );
    BOOST_CHECK( bit_array[3] == 1 );
    BOOST_CHECK( bit_array[4] == 1 );
}

BOOST_AUTO_TEST_CASE( canCreateOneDimensionalBranchAndBoundVariable )
{
    // 14 elements -- or 16 bits (pads to nearest 8 bits)
    haplo::BnbVariable<1> one_dim_variable(14);
    
    // Set element (bit) 2 to 1
    one_dim_variable.set(2, 1);
 
    BOOST_CHECK( one_dim_variable(2) == 1 );   
}

BOOST_AUTO_TEST_SUITE_END()
