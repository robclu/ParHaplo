// ----------------------------------------------------------------------------------------------------------
/// @file   small_containers_gpu.cuh
/// @brief  Header file for small containers for the parahaplo library. This is the equivalent to std::bitset,
///         both in spacial complexity ad in time complexity for access and setting, which can be used with
///         GPUs, unlink bitset
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SMALL_CONTAINERS_GPU_CUH
#define PARAHAPLO_SMALL_CONTAINERS_GPU_CUH

#include <nppdefs.h>

namespace haplo {

using byte = uint8_t;

// ----------------------------------------------------------------------------------------------------------
/// @struct     BinaryOnes
/// @brief      Generates a seqeunce of binary ones at compile time, the length depends on the data type
/// @tparam     SizeType    The type used to determine the size -- not implemented at the momemnt
// ----------------------------------------------------------------------------------------------------------
template <typename SizeType>
struct BinaryOnes {
    __host__ __device__ static constexpr size_t num_elements = sizeof(SizeType) * 8;
    __host__ __device__ static constexpr SizeType value = NPP_MAX_8U ^ (1 << (num_elements - 1));
};

// ----------------------------------------------------------------------------------------------------------
/// @class      TinyContainer
/// @brief      Contianer for binary variables which uses 1 bit pre varable and can hold up to 64 binary
///             variables (DON'T use long or unsigned int -- they aren't working fir bit removal
///             so use byte (8 bits), short (16 bits), int (32 bits), unsigned long (64 bits)
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
    __host__ __device__  
    static constexpr size_t num_elements = sizeof(SizeType) * 8;    
public: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    __host__ __device__
    TinyContainer() : _bits(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bit at position i
    /// @param[in]  i   The index of the bit to get
    /// @return     The value of the bit at position i
    // ------------------------------------------------------------------------------------------------------
    __host__ __device__
    inline byte get(const byte i) const { return (_bits >> (num_elements - 1 - i)) & 0x01; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bit at position i
    /// @param[in]  i       The position of the bit to set
    /// @param[in]  value   The value to set the bit to (0 or 1)
    // ------------------------------------------------------------------------------------------------------
    __host__ __device__
    inline void set(const byte i, const byte value) 
    { 
        // If the value is valid and values != value, then we must swap the values
        if (value <= 1 && (((_bits >> (num_elements - 1 - i)) & 0x01) ^ (value & 0x01))) 
            _bits ^= (0x01 << (num_elements - 1 - i));
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes a bit from the container, and shifts all other bits left, so essentially it's like
    ///             removing and element from the container.
    /// @param[in]  i   The index of the bit to remove
    // ------------------------------------------------------------------------------------------------------
    __host__ __device__
    inline void remove_bit(const byte i) 
    {
        // All ones -- works for signed and unsigned types
        static constexpr SizeType ones = BinaryOnes<SizeType>::value;
        _bits = (_bits & (ones - ((1 << (num_elements - i)) - 1)))     // Get bits to the left of the bit to remove
                ^ (( _bits & ((1 << (num_elements - i - 1)) - 1)) << 1);         // Get bits below the bits to move
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Shifts all bits left by n 
    /// @param[in]  n   THe amount to shift the bits left, default to 1
    // ------------------------------------------------------------------------------------------------------
    __host__ __device__
    inline void shift_left(const size_t n = 1) { _bits <<= n; }
};

   
}
#endif
