// ----------------------------------------------------------------------------------------------------------
/// @file   cuda_link_tests.cpp
/// @brief  Test suite for parahaplo link container test for cuda
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE CudaLinkTests
#endif
#include <boost/test/unit_test.hpp>

#include <cuda.h>
#include "../haplo/link_v2.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

BOOST_AUTO_TEST_SUITE( LinkGpuSuite )

BOOST_AUTO_TEST_CASE( canCreateLinkWithCudaHostVector )
{
    thrust::host_vector<haplo::LinkGpu> link_vector(10);
    
    link_vector[0].homo_weight() = 1;
    link_vector[0].hetro_weight() = 2;
    
    BOOST_CHECK( link_vector[0].homo_weight()  == 1 );
    BOOST_CHECK( link_vector[0].hetro_weight() == 2 );
}

BOOST_AUTO_TEST_CASE( canCreateLinkWithCudaDeviceVector )
{
    thrust::host_vector<haplo::LinkGpu> link_vector(10);
    
    link_vector[0].homo_weight() = 1;
    link_vector[0].hetro_weight() = 2;

    // Check that the links can be assingnem with a device vector    
    thrust::device_vector<haplo::LinkGpu> dlink_vector = link_vector;
}

BOOST_AUTO_TEST_SUITE_END()
