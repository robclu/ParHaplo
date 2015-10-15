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

static constexpr const char* input_four  = "input_files/input_four.txt";

BOOST_AUTO_TEST_SUITE( TreeGpuSuite )

BOOST_AUTO_TEST_CASE( canCreateTree )
{
    using block_type    = haplo::Block<1658, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_four);
    subblock_type   sub_block(block, 1);

    haplo::TreeGpu tree(sub_block.data(), sub_block.read_info(), sub_block.snp_info(), 100);
}

BOOST_AUTO_TEST_SUITE_END()
