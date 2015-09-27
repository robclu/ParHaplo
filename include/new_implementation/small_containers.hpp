// ----------------------------------------------------------------------------------------------------------
/// @file   small_containers.hpp
/// @brief  Header file for small containers for the parahaplo library. This is the equivalent to std::bitset,
///         both in spacial complexity ad in time complexity for access and setting, which can be used with
///         GPUs, unlink bitset
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SMALL_CONTAINERS_HPP
#define PARAHAPLO_SMALL_CONTAINERS_HPP

#include <iostream>
#include <bitset>
#include <limits>

namespace haplo {

using byte = uint8_t;


// ----------------------------------------------------------------------------------------------------------
/// @class      TinyContainer
/// @brief      Contianer for binary variables which uses 1 bit pre varable and can hold up to 64 binary
///             variables 
/// @tparam     SizeType        The type used for the size of the container -- byte for 8 bits, short for 16 
///             bits, int for 32, long for 64 bits etc...
/// @tparam     BitsPerElement  The number of bits per element. The default is 1 bit, but 2 can be used
// ----------------------------------------------------------------------------------------------------------
template <typename SizeType, byte BitsPerElement>
class TinyContainer;

// Specialize for 1 bit per element
template <typename SizeType>
class TinyContainer<SizeType, 1> {
private:
    SizeType                _bits;                                  //!< The bits for the container  
public:
    static constexpr size_t num_elements = sizeof(SizeType) * 8;    //!< The number of bits in the container
public: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    TinyContainer() : _bits(0) {};
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bit at position i
    /// @param[in]  i   The index of the bit to get
    /// @return     The value of the bit at position i
    // ------------------------------------------------------------------------------------------------------
    inline byte get(const byte i) const { return (_bits >> i) & 0x01; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bit at position i
    /// @param[in]  i       The position of the bit to set
    /// @param[in]  value   The value to set the bit to (0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const byte i, const byte value) 
    { 
        // If the value is valid and values != value, then we must swap the values
        if (value <= 1 && ((_bits >> i & 0x01) ^ (value & 0x01))) _bits ^= (0x01 << i);
    }
    
    inline void remove_bit(const byte i) 
    {
        //SizeType ones = (1 << (sizeof(SizeType) * 8)) - 1;
        //SizeType oth = ((1 << (sizeof(SizeType) * 8 - 1)) - 1);
        std::bitset<sizeof(SizeType)*8> a(std::numeric_limits<SizeType>::max() ^ (1 << (sizeof(SizeType) * 8 - 1)));
        //std::bitset<sizeof(SizeType)*8> a(((1 << (sizeof(SizeType) * 8)) - 1) - ((1 << (i + 1)) - 1));
        //std::bitset<sizeof(SizeType)*8> a((_bits & ((1 << i) - 1)) << 1);
        //std::bitset<sizeof(SizeType)*8> a(oth);
        std::cout << a << "\n" << /*b <<*/ "\n";
            
        _bits = (_bits & (0xFF - ((1 << (i + 1)) - 1)))        // Get high bits
                ^ ( _bits & ((1 << i) - 1) << 1);
    }
    
    void print()
    {
        std::bitset<sizeof(SizeType) * 8> x(_bits);
        std::cout << x << "\n";
    }
};

// Specialization for 2 bits per container
template <typename SizeType>
class TinyContainer<SizeType, 2> {
private:
    SizeType                _bits;                                  //!< The bits for the container
public:
    static constexpr size_t num_elements = sizeof(SizeType) * 4;    //!< The number of bits in the container
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    TinyContainer() : _bits(0) {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bits at position i
    /// @param[in]  i   The index of the bits to get
    /// @return     The value of the bits at position i
    // ------------------------------------------------------------------------------------------------------
    inline byte get(const byte i) const { return (_bits >> (i * 2)) & 0x03; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bits at position i
    /// @param[in]  i       The position of the bits to set
    /// @param[in]  value   The value to set the bits to (0 | 1 | 2 | 3)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const byte i, const byte value)
    {
        // Check if the lower (right) bit has the correct value, otherwise set it
        if (value <= 3 && ((_bits >> (i * 2) & 0x01) ^ (value & 0x01))) _bits ^= (0x01 << (i * 2));
        // Check if the higher (left) bit has the correct value, otherwise set it
        if (value <= 3 && ((_bits >> (i * 2 + 1) & 0x01) ^ (value >> 1 & 0x01))) _bits ^= (0x01 << (i * 2 + 1));
    }
};

// ----------------------------------------------------------------------------------------------------------
/// @class  BinaryContainer 
/// @brief  Conatiner which can hold N binary variables, which is optimized for space and performance.
/// @tparam NumElements     The number of binary elements the container can hold
/// @tparam BitsPerElement  The number of bits per element, can be 1 or 2 -- default to 1
// ----------------------------------------------------------------------------------------------------------
template <size_t NumElements, byte BitsPerElement = 1>
class BinaryContainer {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using internal_container = TinyContainer<byte, BitsPerElement>;
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t elements_per_bin    = internal_container::num_elements;
    static constexpr size_t bins                = NumElements / elements_per_bin;;
private:
    internal_container      _data[bins + 1];    //!< Array of bit containers
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    BinaryContainer() {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the binary container
    /// @param[in]  i   The index of the binary in the container to get
    // ------------------------------------------------------------------------------------------------------
    inline byte get(const size_t i) const { return _data[i / elements_per_bin].get(i % elements_per_bin); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const size_t i, const byte value) { _data[i / elements_per_bin].set(i % elements_per_bin, 
                                                                                        value               ); }
};


}           // End namespace haplo

#endif      // PARAHAPLO_VARIABLES_HPP

