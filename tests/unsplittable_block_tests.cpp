// ----------------------------------------------------------------------------------------------------------
/// @file   unsplittable_block_tests.cpp
/// @brief  Test suite for parahaplo unsplittable block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/unsplittable_block_cpu.hpp"

static constexpr const char* input_1 = "input_files/input_unfiltered_1.txt";
static constexpr const char* input_2 = "input_files/input_unfiltered_2.txt";
static constexpr const char* input_3 = "input_files/input_singleton_rows.txt";

BOOST_AUTO_TEST_SUITE( UnsplittableBlockSuite )
   
// NOTE: This test doesn't actually test anything, but the outputs of 
//       the test should show that an out of range exception was thrown
BOOST_AUTO_TEST_CASE( errorIsThrownForOutOfRangeUnsplittableBlock  )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> unsplittable_block(block, 7);
}

BOOST_AUTO_TEST_CASE( canCreateUnsplittableBlockCorrectlyAndGetData )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> unsplittable_block(block, 0);
}



BOOST_AUTO_TEST_SUITE_END()
