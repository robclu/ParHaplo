// ----------------------------------------------------------------------------------------------------------
/// @file   small_containers.hpp
/// @brief  Header file for small containers for the parahaplo library. This is the equivalent to std::bitset,
///         both in spacial complexity ad in time conplexity for access and setting
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SMALL_CONTAINERS_HPP
#define PARAHAPLO_SMALL_CONTAINERS_HPP

namespace haplo {

using byte = uint8_t;

template <typename SizeType>
class TinyContainer {
public:
    SizeType _bits;
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
    inline SizeType get(const byte i) const { return (_bits >> i) & 0x01; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bit at position i
    /// @param[in]  i       THe position of the bit to set
    /// @param[in]  value   The value to set the bit to (0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const byte i, const byte value) 
    { 
        // If the value is valid and values != value, then we must swap the values
        if (value <= 1 && !((_bits >> i & 0x01) && value)) _bits ^= (0x01 << i);
    }
};

// ----------------------------------------------------------------------------------------------------------
/// @class  BinaryContainer 
/// @brief  Conatiner which can hold N binary variables, which is optimized for space and performance.
/// @tparam NumElements     The number of binary elements the container can hold
// ----------------------------------------------------------------------------------------------------------
template <size_t NumElements>
class BinaryContainer {
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    BinaryContainer() {}
    
    // Define the bin offset
    static constexpr size_t bins = NumElements / 8;
private:
    TinyContainer<byte> _data[bins + 1];
    
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the binary container
    /// @param[in]  i   The index of the binary in the container to get
    // ------------------------------------------------------------------------------------------------------
    inline byte get(const size_t i) const { return _data[i / 8].get(i % 8); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const size_t i, const byte value) { _data[i / 8].set(i % 8, value); }
};

}           // End namespace haplo

#endif      // PARAHAPLO_VARIABLES_HPP

