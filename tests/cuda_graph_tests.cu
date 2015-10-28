// ----------------------------------------------------------------------------------------------------------
/// @file   cuda_graph_tests.cu
/// @brief  Test suite for parahaplo GPU graph search tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE CudaGraphTests
#endif
#include <boost/test/unit_test.hpp>
#include <chrono>

#include "../haplo/subblock_cpu.hpp"
#include "../haplo/graph_gpu.h"
#include "cuda_runtime.h"

using namespace std::chrono;

static constexpr const char* test_input    = "output_files/output_simulated_265.txt";
static constexpr const char* test_input1   = "new_outputs/geraci/100_10_0.1_0.4/output_5298.txt";

BOOST_AUTO_TEST_SUITE( GraphGpuSuite )

BOOST_AUTO_TEST_CASE( canCreateGraph )
{
    using block_type    = haplo::Block<5298, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    using graph_type    = haplo::Graph<subblock_type, haplo::devices::gpu>;

    // System timer
    high_resolution_clock::time_point start = high_resolution_clock::now(); 

    block_type block(test_input1);
    
    high_resolution_clock::time_point end = high_resolution_clock::now();
    duration<double> sort_time            = duration_cast<duration<double>>(end - start);
    std::cout << sort_time.count() << "\n";
    
    subblock_type sub_block(block, 1);

    end         = high_resolution_clock::now();
    sort_time   = duration_cast<duration<double>>(end - start);
    std::cout << sort_time.count() << "\n";
    
    // For now just use the first device 
    // Later use device manager
    size_t device_index = 0;
    
    graph_type graph(sub_block, device_index);
    
    // Search the tree for the haplotype 
    graph.search();
    
    end         = high_resolution_clock::now();
    sort_time   = duration_cast<duration<double>>(end - start);
    std::cout << sort_time.count() << "\n";
    
    // Put the haplotypes back into the block
    block.merge_haplotype(sub_block);
    
    end         = high_resolution_clock::now();
    sort_time   = duration_cast<duration<double>>(end - start);
    std::cout << sort_time.count() << "\n";
    
    graph.print_mec();
    block.print_haplotypes();
}

BOOST_AUTO_TEST_SUITE_END()
