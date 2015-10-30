// ----------------------------------------------------------------------------------------------------------
/// @file   subblock_tests.cpp
/// @brief  Test suite for parahaplo sub-block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SubBlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/subblock_cpu.hpp"
#include <chrono>

using namespace std::chrono;

static constexpr const char* input_zero   = "input_files/input_zero.txt";
static constexpr const char* input_one    = "input_files/input_one.txt";
static constexpr const char* input_two    = "input_files/input_two.txt";
static constexpr const char* input_three  = "input_files/input_three.txt";
static constexpr const char* input_four  = "input_files/input_four.txt";

BOOST_AUTO_TEST_SUITE( SubBlockSuite )

// NOTE: This test doesn't actually test anything, but the outputs of 
//       the test should show that an out of range exception was thrown
BOOST_AUTO_TEST_CASE( errorIsThrownForOutOfRangeSubBlock  )
{
    // Define a block for a with 4 CPU cores
    using block_type = haplo::Block<28, 2, 2>;
    
    // First create the block
    block_type block(input_zero);
    
    // Create a sub-block from the block out of range -- just to illustrate out of range error
    haplo::SubBlock<block_type, 2, 2, haplo::devices::cpu> sub_block(block, 7);
}

BOOST_AUTO_TEST_CASE( canCreateSubBlockCorrectlyAndGetData1 )
{
    // Define a block for a with 4 CPU cores
    using block_type = haplo::Block<28, 2, 2>;
    
    // First create the block
    block_type block(input_zero);

    haplo::SubBlock<block_type, 2, 2, haplo::devices::cpu> subblock(block, 0);

    BOOST_CHECK( subblock(0, 0) == 0 );
    BOOST_CHECK( subblock(0, 1) == 0 );
    BOOST_CHECK( subblock(0, 2) == 3 );
    BOOST_CHECK( subblock(1, 0) == 1 );
    BOOST_CHECK( subblock(1, 1) == 0 );
    BOOST_CHECK( subblock(1, 2) == 0 );
    BOOST_CHECK( subblock(2, 0) == 1 );
    BOOST_CHECK( subblock(2, 1) == 1 );
    BOOST_CHECK( subblock(2, 2) == 3 );
}

BOOST_AUTO_TEST_CASE( canRemoveMonotoneColumns )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2);
              
    BOOST_CHECK( sub_block.base_start_row() == 4 ); 
    BOOST_CHECK( sub_block(0, 0)  == 3 );
    BOOST_CHECK( sub_block(0, 1)  == 3 );
    BOOST_CHECK( sub_block(0, 2)  == 0 );
    BOOST_CHECK( sub_block(0, 3)  == 0 );
    BOOST_CHECK( sub_block(1, 0)  == 3 );
    BOOST_CHECK( sub_block(1, 1)  == 3 );
    BOOST_CHECK( sub_block(1, 2)  == 1 );
    BOOST_CHECK( sub_block(1, 3)  == 1 );
    BOOST_CHECK( sub_block(2, 0)  == 0 );
    BOOST_CHECK( sub_block(2, 1)  == 1 );
    BOOST_CHECK( sub_block(2, 2)  == 2 );
    BOOST_CHECK( sub_block(2, 3)  == 1 );
    BOOST_CHECK( sub_block(3, 0)  == 3 );
    BOOST_CHECK( sub_block(3, 1)  == 0 );
    BOOST_CHECK( sub_block(3, 2)  == 1 );
    BOOST_CHECK( sub_block(3, 3)  == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
