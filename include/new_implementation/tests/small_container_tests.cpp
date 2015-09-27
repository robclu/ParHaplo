// ----------------------------------------------------------------------------------------------------------
/// @file   small_container_tests.cpp
/// @brief  Test suite for parahaplo small container tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SmallContainerTests
#endif
#include <boost/test/unit_test.hpp>

#include "../small_containers.hpp"


BOOST_AUTO_TEST_SUITE( BinaryContainerSuite )
    
BOOST_AUTO_TEST_CASE( canCreateBinaryContainerWith1BitPerElement ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<12> elements;
    
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

BOOST_AUTO_TEST_CASE( canCreateBinaryContainerWith2BitsPerElement ) 
{
    // Define the conatiner to use 2 bits per element and 
    // 18 elements so 36 bits = 5 bytes
    haplo::BinaryContainer<18, 2> elements;

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

BOOST_AUTO_TEST_CASE( canRemoveBitsFrom1BitBinaryContainer ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::TinyContainer<haplo::byte, 1> bits;
    
    // Set some elements
    bits.set(1, 1);
    bits.set(2, 1);
    bits.set(7, 1);
    
    bits.remove_bit(7);

    BOOST_CHECK( bits.get(7) == 0 );
    BOOST_CHECK( bits.get(1) == 0 );
    BOOST_CHECK( bits.get(2) == 1 );
    BOOST_CHECK( bits.get(3) == 1 );
}

BOOST_AUTO_TEST_CASE( canRemoveBitsFrom2BitBinaryContainer ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::TinyContainer<haplo::byte, 2> bits;
    
    // Set some elements
    bits.set(0, 3);
    bits.set(1, 1);
    bits.set(2, 2);
    
    bits.remove_bit(2);
    
    bits.shift_left(1);
    
    BOOST_CHECK( bits.get(0) == 0 );
    BOOST_CHECK( bits.get(1) == 0 );
    BOOST_CHECK( bits.get(2) == 3 );
    BOOST_CHECK( bits.get(3) == 2 );
}

BOOST_AUTO_TEST_CASE( canRemoveElementsOf1BitContainer ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<12> elements;
    
    // Set some elements
    elements.set(1, 1);
    elements.set(2, 1);
    elements.set(6, 1);
    elements.set(7, 1);
    elements.set(9, 1);    
    elements.set(3, 0);
    
    elements.remove_element(6);
    elements.remove_element(7);

    BOOST_CHECK( elements.get(1) == 0 );
    BOOST_CHECK( elements.get(2) == 0 );
    BOOST_CHECK( elements.get(3) == 1 );
    BOOST_CHECK( elements.get(4) == 1 );
    BOOST_CHECK( elements.get(5) == 0 );
    BOOST_CHECK( elements.get(6) == 0 );
    BOOST_CHECK( elements.get(7) == 0 );
    BOOST_CHECK( elements.get(9) == 1 );
    BOOST_CHECK( elements.size() == 10 );
}

BOOST_AUTO_TEST_CASE( canRemoveElementsOf12BitContainer ) 
{
    // Define the binary container to use 2 bits per element
    haplo::BinaryContainer<14, 2> elements;
    
    // Set some elements
    elements.set(1, 1);
    elements.set(2, 2);
    elements.set(6, 0);
    elements.set(7, 3);
    elements.set(9, 2);    
    elements.set(3, 0);

    elements.print();

    std::cout << "\n";    
    elements.remove_element(6);
    
    elements.print(); std::cout << "\n";
    
    elements.remove_element(2);

    BOOST_CHECK( elements.get(1) == 0 );
    BOOST_CHECK( elements.get(2) == 0 );
    BOOST_CHECK( elements.get(3) == 1 );
    BOOST_CHECK( elements.get(4) == 0 );
    BOOST_CHECK( elements.get(5) == 0 );
    BOOST_CHECK( elements.get(6) == 0 );
    BOOST_CHECK( elements.get(7) == 3 );
    BOOST_CHECK( elements.get(9) == 2 );
    BOOST_CHECK( elements.size() == 12 );
    
    elements.print();
}

BOOST_AUTO_TEST_SUITE_END()
