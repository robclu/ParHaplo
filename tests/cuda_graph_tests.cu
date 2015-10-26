// ----------------------------------------------------------------------------------------------------------
/// @file   cuda_graph_tests.cu
/// @brief  Test suite for parahaplo GPU graph search tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE CudaGraphTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/subblock_cpu.hpp"
#include "../haplo/graph_gpu.h"
#include "cuda_runtime.h"

using namespace std::chrono;

static constexpr const char* test_input  = "graph_inputs/test_input.txt";
static constexpr const char* input_zero   = "input_files/input_zero.txt";

BOOST_AUTO_TEST_SUITE( GraphGpuSuite )

BOOST_AUTO_TEST_CASE( canCreateGraph )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    using graph_type    = haplo::Graph<subblock_type, haplo::devices::gpu>;
   
    block_type block(input_zero);
    
    std::cout << "NUM_SUB_BLOCKS " <<  block.num_subblocks() << "\n";

    subblock_type sub_block(block, 0);
    sub_block.print();

    // For now just use the first device 
    // Later use device manager
    size_t device_index = 0;
    
    graph_type graph(sub_block, device_index);
    
    // Search the tree for the haplotype 
    graph.search();
}

BOOST_AUTO_TEST_SUITE_END()
