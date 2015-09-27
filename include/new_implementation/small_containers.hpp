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
/// @struct     BinaryOnes
/// @brief      Generates a seqeunce of binary ones at compile time, the length depends on the data type
/// @tparam     SizeType    The type used to determine the size
// ----------------------------------------------------------------------------------------------------------
template <typename SizeType>
struct BinaryOnes {
    static constexpr size_t num_elements = sizeof(SizeType) * 8;
    static constexpr SizeType value = std::numeric_limits<SizeType>::max() ^ (1 << (num_elements - 1));
};

// specialization for 8 bit unsigned
template <>
struct BinaryOnes<byte> {
    static constexpr byte value = std::numeric_limits<byte>::max();
};

// ----------------------------------------------------------------------------------------------------------
/// @class      TinyContainer
/// @brief      Contianer for binary variables which uses 1 bit pre varable and can hold up to 64 binary
///             variables (DON'T use long or unsigned int -- they aren't working fir bit removal
///             so use byte (8 bits), short (16 bits), int (32 bits), unsigned lon (64 bits)
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

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes a bit from the container, and shifts all other bits left, so essentially it's like
    ///             removing and element from the container.
    /// @param[in]  i   The index of the bit to remove
    // ------------------------------------------------------------------------------------------------------
    inline void remove_bit(const byte i) 
    {
        // All ones -- works for signed and unsigned types
        static constexpr SizeType ones = BinaryOnes<SizeType>::value;
        _bits = (_bits & (ones - ((1 << (i + 1)) - 1)))     // Get bits to the left of the bit to remove
                ^ (( _bits & ((1 << i) - 1)) << 1);         // Get bits below the bits to move
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Shifts all bits left by n 
    /// @param[in]  n   THe amount to shift the bits left, default to 1
    // ------------------------------------------------------------------------------------------------------
    inline void shift_left(const size_t n = 1) { _bits <<= n; }
    
    void print()
    {
        std::bitset<sizeof(SizeType) * 8> x(_bits);
        std::cout << x;
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

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes a bit from the container, and shifts all other bits left, so essentially it's like
    ///             removing and element from the container.
    /// @param[in]  i   The index of the bit to remove
    // ------------------------------------------------------------------------------------------------------
    inline void remove_bit(const byte i) 
    {
        // All ones -- works for signed and unsigned types
        static constexpr SizeType ones = BinaryOnes<SizeType>::value;
        _bits = (_bits & (ones - ((1 << (2 * i + 2)) - 1)))     // Get bits to the left of the bit to remove
                ^ (( _bits & ((1 << (2 * i)) - 1)) << 2);        // Get bits below the bits to move
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Shifts all bits left by n 
    /// @param[in]  n   THe amount (number of positions) to shift the bits left, default to 1
    // ------------------------------------------------------------------------------------------------------
    inline void shift_left(const size_t n = 1) { _bits <<= (2 * n); }
    
    void print()
    {
        std::bitset<num_elements * 2> x(_bits);
        std::cout << x;
    }
};

// ----------------------------------------------------------------------------------------------------------
/// @class  BinaryContainer 
/// @brief  Container which can hold N binary variables, which is optimized for space and performance. The
///         bits are stored little endian
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
    size_t                  _num_elements;      //!< Number of elements in the container
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    BinaryContainer() : _num_elements(NumElements) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the binary container
    /// @param[in]  i   The index of the binary in the container to get
    // ------------------------------------------------------------------------------------------------------
    inline byte get(const size_t i) const 
    { 
        return _data[bins - (i / elements_per_bin)].get(i % elements_per_bin); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const size_t i, const byte value) 
    { 
        _data[bins - (i / elements_per_bin)].set(i % elements_per_bin, value); 
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The size (number of elements) in the container 
    /// @return     The number of elements in the container 
    // ------------------------------------------------------------------------------------------------------
    inline size_t size() const { return _num_elements; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes an element from the container
    /// @param[in]  i   The index of the element to remove
    // ------------------------------------------------------------------------------------------------------
    inline void remove_element(const size_t i) 
    {
        // Decrase the number of elements
        --_num_elements;
        
        // First remove the bit which is the element
        _data[bins - (i / elements_per_bin)].remove_bit(i % elements_per_bin);

        // Now for each of the bins to the right, set the LSB of the container to the MSB of the
        // container to the right, then shift all bits of the container to the right, to the left
        for (size_t bin = bins - (i / elements_per_bin); bin < bins; ++bin) {
            _data[bin].set(0, _data[bin + 1].get(elements_per_bin - 1));
            _data[bin + 1].shift_left(1);
        }
    }
    
    void print() 
    {
        for (int i = 0; i < bins + 1; ++i)
            _data[i].print();
    }
};


}           // End namespace haplo

#endif      // PARAHAPLO_VARIABLES_HPP

