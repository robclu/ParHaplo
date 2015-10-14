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
#include <cuda_runtime.h>
#include <thrust/device_vector.h>

    struct TestB {
        int x,y,z;
        TestB() : x(1), y(2), z(3) {}
        
    };
    
    struct TestA {
        int x, y;
        
        TestA(TestB b) : x(b.x), y(b.y) {}
    };

BOOST_AUTO_TEST_SUITE( BinaryArraySuite )
 
// --------------------------------- TINY CONTAINER TESTS ---------------------------------------------------

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
    
    
    thrust::host_vector<TestB> tb(20);
    
    thrust::device_vector<TestA> ta = tb;
}

// ------------------------------------------ VECTOR TESTS --------------------------------------------------

BOOST_AUTO_TEST_CASE( canCreateBinaryVectorWith1BitPerElement ) 
{
    // BinaryVector uses 1 bit per element
    haplo::BinaryVector<1> elements(12);
    
    // Set some elements
    elements.set(1, 1);
    elements.set(2, 1);
    elements.set(7, 1);
    elements.set(9, 1);    
    elements.set(3, 0);
    
    BOOST_CHECK( elements.get(1) == 1 );
    BOOST_CHECK( elements.get(2) == 1 );
    BOOST_CHECK( elements.get(3) == 0 );
    BOOST_CHECK( elements.get(7) == 1 );
    BOOST_CHECK( elements.get(8) == 0 );
    BOOST_CHECK( elements.get(9) == 1 );
}

BOOST_AUTO_TEST_CASE( canCreateBinaryVectorWith2BitsPerElement ) 
{
    // Define the conatiner to use 2 bits per element and 
    haplo::BinaryVector<2> elements(18);

    // Set some elements
    elements.set(1 , 1);
    elements.set(3 , 2);
    elements.set(7 , 1);
    elements.set(14, 3);
    elements.set(8 , 0);
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 0 );
    BOOST_CHECK( elements.get(3)  == 2 );
    BOOST_CHECK( elements.get(7)  == 1 );
    BOOST_CHECK( elements.get(8)  == 0 );
    BOOST_CHECK( elements.get(14) == 3 );
}

BOOST_AUTO_TEST_CASE( canRemoveElementsOf1BitBinaryVector ) 
{
    haplo::BinaryVector<1> elements(12);
    
    // Set some elements
    elements.set(1, 1);
    elements.set(2, 1);
    elements.set(6, 1);
    elements.set(7, 1);
    elements.set(9, 1);    
    elements.set(3, 0);

    elements.remove_element(6);   
    elements.remove_element(7);
    
    BOOST_CHECK( elements.get(1) == 1 );
    BOOST_CHECK( elements.get(2) == 1 );
    BOOST_CHECK( elements.get(3) == 0 );
    BOOST_CHECK( elements.get(4) == 0 );
    BOOST_CHECK( elements.get(5) == 0 );
    BOOST_CHECK( elements.get(6) == 1 );
    BOOST_CHECK( elements.get(7) == 1 );
    BOOST_CHECK( elements.get(9) == 0 );
    BOOST_CHECK( elements.size() == 10 );
}

BOOST_AUTO_TEST_CASE( canRemoveElementsOf2BitBinaryVector ) 
{
    haplo::BinaryVector<2> elements(14);
    
    // Set some elements
    elements.set(1 , 1);
    elements.set(2 , 2);
    elements.set(3 , 0);
    elements.set(6 , 0);
    elements.set(7 , 3);
    elements.set(9 , 2);    
    elements.set(11, 3);

    // Elements looks like : 00011000000000110010001100000000
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 2 );
    BOOST_CHECK( elements.get(3)  == 0 );
    BOOST_CHECK( elements.get(6)  == 0 );
    BOOST_CHECK( elements.get(7)  == 3 );
    BOOST_CHECK( elements.get(9)  == 2 );
    BOOST_CHECK( elements.get(11) == 3 );
    BOOST_CHECK( elements.get(12) == 0 );
    BOOST_CHECK( elements.size() == 14 );
    
    elements.remove_element(6);

    // Now elements looks like : 0001100000001100100011000000
    
    elements.remove_element(2);

    // Now elements looks like : 0001000000110010001100000000

    BOOST_CHECK( elements.get(0) == 0 );
    BOOST_CHECK( elements.get(1) == 1 );
    BOOST_CHECK( elements.get(2) == 0 );
    BOOST_CHECK( elements.get(3) == 0 );
    BOOST_CHECK( elements.get(4) == 0 );
    BOOST_CHECK( elements.get(5) == 3 );
    BOOST_CHECK( elements.get(6) == 0 );
    BOOST_CHECK( elements.get(7) == 2 );
    BOOST_CHECK( elements.get(8) == 0 );
    BOOST_CHECK( elements.get(9) == 3 );
    BOOST_CHECK( elements.size() == 12 );
}

BOOST_AUTO_TEST_CASE( canResizeBianryVectorFromZero ) 
{
    haplo::BinaryVector<2> elements;
    
    elements.resize(14);
    
    // Set some elements
    elements.set(1 , 1);
    elements.set(2 , 2);
    elements.set(3 , 0);
    elements.set(6 , 0);
    elements.set(7 , 3);
    elements.set(9 , 2);    
    elements.set(11, 3);

    // Elements looks like : 00011000000000110010001100000000
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 2 );
    BOOST_CHECK( elements.get(3)  == 0 );
    BOOST_CHECK( elements.get(6)  == 0 );
    BOOST_CHECK( elements.get(7)  == 3 );
    BOOST_CHECK( elements.get(9)  == 2 );
    BOOST_CHECK( elements.get(11) == 3 );
    BOOST_CHECK( elements.get(12) == 0 );
    BOOST_CHECK( elements.size() == 14 );
}

BOOST_AUTO_TEST_CASE( canResizeBinaryVectorUpwards )
{
    haplo::BinaryVector<2> elements;
    
    elements.resize(14);
    
    // Set some elements
    elements.set(1 , 1);
    elements.set(2 , 2);
    elements.set(3 , 0);
    elements.set(6 , 0);
    elements.set(7 , 3);
    elements.set(9 , 2);    
    elements.set(11, 3);
    
    elements.resize(18);
    
    elements.set(15, 1);
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 2 );
    BOOST_CHECK( elements.get(3)  == 0 );
    BOOST_CHECK( elements.get(6)  == 0 );
    BOOST_CHECK( elements.get(7)  == 3 );
    BOOST_CHECK( elements.get(9)  == 2 );
    BOOST_CHECK( elements.get(11) == 3 );
    BOOST_CHECK( elements.get(12) == 0 );
    BOOST_CHECK( elements.get(14) == 0 );
    BOOST_CHECK( elements.get(15) == 1 );
    BOOST_CHECK( elements.size() == 18 );
}

BOOST_AUTO_TEST_CASE( canResizeBinaryVectorDownwards ) 
{
    haplo::BinaryVector<2> elements;
    
    elements.resize(14);
    
    // Set some elements
    elements.set(1 , 1);
    elements.set(2 , 2);
    elements.set(3 , 0);
    elements.set(6 , 0);
    elements.set(7 , 3);
    elements.set(9 , 2);    
    elements.set(11, 3);
    
    elements.resize(10);
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 2 );
    BOOST_CHECK( elements.get(3)  == 0 );
    BOOST_CHECK( elements.get(6)  == 0 );
    BOOST_CHECK( elements.get(7)  == 3 );
    BOOST_CHECK( elements.get(9)  == 2 );
    BOOST_CHECK( elements.size() == 10 );
}

BOOST_AUTO_TEST_CASE( canPushBackToVector )
{
    haplo::BinaryVector<2> elements(15);
    
    elements.push_back(1);
    
    BOOST_CHECK( elements.size()    == 16 );
    BOOST_CHECK( elements.get(15)   == 1  );

    elements.push_back(2);    
    elements.push_back(3);    
    
    BOOST_CHECK( elements.size()    == 18 );
    BOOST_CHECK( elements.get(16)   == 2  );
    BOOST_CHECK( elements.get(17)   == 3  );
    
    // Check that the vector did the over-allocation
    BOOST_CHECK( elements.get(18)   == 0  );
}

BOOST_AUTO_TEST_SUITE_END()
