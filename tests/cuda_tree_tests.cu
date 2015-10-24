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
#include "../haplo/tree_gpu.h"
#include "cuda_runtime.h"

using namespace std::chrono;

static constexpr const char* input_100_10  = "new_outputs/geraci/100_10_0.1_0.4/output_5298.txt";

BOOST_AUTO_TEST_SUITE( TreeGpuSuite )

BOOST_AUTO_TEST_CASE( canCreateTree )
{
    using block_type    = haplo::Block<5298, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    using tree_type     = haplo::Tree<subblock_type, haplo::devices::gpu>;
    
    block_type      block(input_100_10);
    subblock_type   sub_block(block, 1);
    sub_block.print();

    // For now just use the first device 
    // Later use device manager
    size_t device_index = 0;
    
    tree_type tree(sub_block, device_index);
    
    // Search the tree for the haplotype 
    tree.search();
}

BOOST_AUTO_TEST_SUITE_END()
