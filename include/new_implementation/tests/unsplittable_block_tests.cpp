// ----------------------------------------------------------------------------------------------------------
/// @file   unsplittable_block_tests.cpp
/// @brief  Test suite for parahaplo unsplittable block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../unsplittable_block_cpu.hpp"

constexpr char* input_1 = "input_unfiltered_1.txt";
constexpr char* input_2 = "input_unfiltered_2.txt";
constexpr char* input_3 = "input_singleton_rows.txt";

BOOST_AUTO_TEST_SUITE( UnsplittableBlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateAnUnsplittableBlock )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> unsplittable_block(block);
}



BOOST_AUTO_TEST_SUITE_END()
