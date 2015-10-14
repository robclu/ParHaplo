// ----------------------------------------------------------------------------------------------------------
/// @file   snp_info_tests.cu
/// @brief  Test suite for parahaplo SnpInfo classes
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE EvaluatorTests
#endif
#include <boost/test/unit_test.hpp>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "../haplo/snp_info_gpu.h"

BOOST_AUTO_TEST_SUITE( SnpSuite )

// Checking for ocmpilation    
BOOST_AUTO_TEST_CASE( canCreateHostAndDeviceVectorsForSnps )
{
    thrust::host_vector<haplo::SnpInfo> host_snps(10);
    
    host_snps[0].start_index() = 4;
    
    thrust::device_vector<haplo::SnpInfoGpu> device_snps = host_snps;
    
    thrust::host_vector<haplo::SnpInfoGpu> host_gsnps = device_snps;
    
    BOOST_CHECK( host_gsnps[0].start_index() == 4 );
}

BOOST_AUTO_TEST_SUITE_END()
