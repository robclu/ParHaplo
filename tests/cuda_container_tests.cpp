// ----------------------------------------------------------------------------------------------------------
/// @file   cuda_container_tests.cpp
/// @brief  Test suite for parahaplo small container tests for cuda
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE CudaContainerTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/small_containers_gpu.cuh"

#include <tbb/concurrent_vector.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

BOOST_AUTO_TEST_SUITE( BinaryArraySuite )
 
// ------------------------------------------ ARRAY TESTS ---------------------------------------------------

BOOST_AUTO_TEST_CASE( canRemoveBitsFrom1BitTinyContainer ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::TinyContainer<haplo::byte, 1> bits;
    
    // Set some elements
    bits.set(1, 1);
    bits.set(2, 1);
    bits.set(7, 1);
    bits.set(3, 0);
    
    // Big endian so bits = 01100001

    BOOST_CHECK( bits.get(1) == 1 );
    BOOST_CHECK( bits.get(2) == 1 );
    BOOST_CHECK( bits.get(3) == 0 );
    BOOST_CHECK( bits.get(7) == 1 );
        
    bits.remove_bit(7);             // Removes bit 7 : 01100000
    bits.remove_bit(0);             // Removes bit 0 : 11000000
    bits.remove_bit(1);             // Removes bit 1 : 10000000
    
    BOOST_CHECK( bits.get(0) == 1 );
    BOOST_CHECK( bits.get(1) == 0 );
    BOOST_CHECK( bits.get(2) == 0 );
    BOOST_CHECK( bits.get(3) == 0 );
    BOOST_CHECK( bits.get(4) == 0 );
    BOOST_CHECK( bits.get(5) == 0 );
    BOOST_CHECK( bits.get(6) == 0 );
    BOOST_CHECK( bits.get(7) == 0 );
    
    tbb::concurrent_vector<size_t> tv(10);
    
    for (size_t i = 0; i < 10; i++) tv.push_back(i);
    

    thrust::host_vector<size_t> dv = std::move(tv);
}

BOOST_AUTO_TEST_SUITE_END()
