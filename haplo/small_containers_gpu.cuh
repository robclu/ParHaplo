// ----------------------------------------------------------------------------------------------------------
/// @file   small_containers_gpu.cuh
/// @brief  Header file for small containers for the parahaplo library. This is the equivalent to std::bitset,
///         both in spacial complexity ad in time complexity for access and setting, which can be used with
///         GPUs, unlink bitset
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SMALL_CONTAINERS_GPU_CUH
#define PARAHAPLO_SMALL_CONTAINERS_GPU_CUH

#include <nppdefs.h>
#include <thrust/host_vector.h>

#ifdef __CUDACC__
    #define CUDA_HD __host__ __device__
    #define CUDA_H  __host__
    #define CUDA_D  __device__
    #define SHARED  __shared__
#else
    #define CUDA_HD
    #define CUDA_H
    #define CUDA_D
    #define SHARED
#endif

namespace haplo {

using byte = uint8_t;

// ----------------------------------------------------------------------------------------------------------
/// @struct     BinaryOnes
/// @brief      Generates a seqeunce of binary ones at compile time, the length depends on the data type
/// @tparam     SizeType    The type used to determine the size -- not implemented at the momemnt
// ----------------------------------------------------------------------------------------------------------
template <typename SizeType>
struct BinaryOnes {
    static const size_t num_elements = sizeof(SizeType) * 8;
    static const  SizeType value = NPP_MAX_8U ^ (1 << (num_elements - 1));
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
    static const size_t num_elements = sizeof(SizeType) * 8;    
public: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    TinyContainer() : _bits(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bit at position i
    /// @param[in]  i   The index of the bit to get
    /// @return     The value of the bit at position i
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline byte get(const byte i) const { return (_bits >> (num_elements - 1 - i)) & 0x01; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bit at position i
    /// @param[in]  i       The position of the bit to set
    /// @param[in]  value   The value to set the bit to (0 or 1)
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
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
    CUDA_HD
    inline void remove_bit(const byte i) 
    {
        // All ones -- works for signed and unsigned types
        const SizeType ones = BinaryOnes<SizeType>::value;
        _bits = (_bits & (ones - ((1 << (num_elements - i)) - 1)))     // Get bits to the left of the bit to remove
                ^ (( _bits & ((1 << (num_elements - i - 1)) - 1)) << 1);         // Get bits below the bits to move
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Shifts all bits left by n 
    /// @param[in]  n   THe amount to shift the bits left, default to 1
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void shift_left(const size_t n = 1) { _bits <<= n; }
};

// Specialization for 2 bits per container
template <typename SizeType>
class TinyContainer<SizeType, 2> {
private:
    SizeType                _bits;                                  //!< The bits for the container
public:
    using standard_container = thrust::host_vector<uint8_t>;
    static const size_t num_elements = sizeof(SizeType) * 4;    //!< The number of bits in the container
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    TinyContainer() : _bits(0) {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bits at position i
    /// @param[in]  i   The index of the bits to get
    /// @return     The value of the bits at position i
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline byte get(const byte i) const { return (_bits >> ((num_elements - i - 1) * 2)) & 0x03; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bits at position i
    /// @param[in]  i       The position of the bits to set
    /// @param[in]  value   The value to set the bits to (0 | 1 | 2 | 3)
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void set(const byte i, const byte value)
    {
        // Check if the lower (right) bit has the correct value, otherwise set it
        if (value <= 3 && ((_bits >> ((num_elements - i - 1) * 2) & 0x01) ^ (value & 0x01))) 
            _bits ^= (1 << ((num_elements - i - 1) * 2));
        // Check if the higher (left) bit has the correct value, otherwise set it
        if (value <= 3 && ((_bits >> ((num_elements - i - 1) * 2 + 1) & 0x01) ^ (value >> 1 & 0x01))) 
            _bits ^= (0x01 << ((num_elements - i - 1) * 2 + 1));
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes a bit from the container, and shifts all other bits left, so essentially it's like
    ///             removing and element from the container.
    /// @param[in]  i   The index of the bit to remove
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void remove_bit(const byte i) 
    {
        // All ones -- works for signed and unsigned types
        const SizeType ones = BinaryOnes<SizeType>::value;
        _bits = (_bits & (ones - ((1 << ((num_elements - i) * 2)) - 1)))     // Get bits to the left of the bit to remove
                ^ (( _bits & ((1 << ((num_elements - i - 1) * 2)) - 1)) << 2);        // Get bits below the bits to move
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Shifts all bits left by n 
    /// @param[in]  n   THe amount (number of positions) to shift the bits left, default to 1
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void shift_left(const size_t n = 1) { _bits <<= (2 * n); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts to bool vector
    // ------------------------------------------------------------------------------------------------------
    thrust::host_vector<uint8_t>  to_thrust_vector() const 
    {
        standard_container host_vec;
        for (auto i = 0; i < num_elements; ++i)
            host_vec.push_back(get(i));
        return host_vec;
    }
}; 

// ----------------------------------------------------------------------------------------------------------
/// @class  BinaryVector
/// @brief  Container which can hold N binary variables, which is optimized for space and performance. The
///         bits are stored big endian
/// @tparam BitsPerElement  The number of bits per element, can be 1 or 2 -- default to 1
// ----------------------------------------------------------------------------------------------------------
template <byte BitsPerElement = 1>
class BinaryVector {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using standard_container = thrust::host_vector<uint8_t>;
    using internal_container = TinyContainer<byte, BitsPerElement>;
    using data_container     = thrust::host_vector<internal_container>;
    // ------------------------------------------------------------------------------------------------------
    static const size_t elements_per_bin    = internal_container::num_elements;
private:
    data_container          _data;              //!< Vector of bit containers
    size_t                  _bins;              //!< The number of bins in the vector
    size_t                  _num_elements;      //!< Number of elements in the container
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    explicit BinaryVector() : _data(0), _bins(0), _num_elements(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to take the number of elements for the container 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    explicit BinaryVector(const size_t num_elements) 
    : _data(num_elements / elements_per_bin + 1), _bins(num_elements / elements_per_bin), 
      _num_elements(num_elements)
    {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the binary container
    /// @param[in]  i   The index of the binary in the container to get
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline byte get(const size_t i) const 
    { 
        return _data[i / elements_per_bin].get(i % elements_per_bin); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void set(const size_t i, const byte value) 
    { 
        _data[i / elements_per_bin].set(i % elements_per_bin, value); 
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The size (number of elements) in the container 
    /// @return     The number of elements in the container 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t size() const { return _num_elements; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resizes the data container (adds to the end if bigger, or removes from the end if smaller) 
    /// @param[in]  num_elements    The number of elements in the container after the resize
    // ------------------------------------------------------------------------------------------------------
    CUDA_H
    inline void resize(const size_t num_elements) 
    { 
        size_t total_bins       = num_elements / elements_per_bin + 1;
        size_t current_bins     = _bins;
        
        if (total_bins > current_bins) {
            for (size_t i = current_bins; i < total_bins; ++i) {
                _data.push_back(internal_container());
                ++_bins;
            }
        } else if (total_bins < current_bins) {
            for (size_t i = current_bins; i > total_bins; --i) {
                _data.pop_back();
                --_bins;
            }
        }
        _num_elements = num_elements;
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes an element from the container
    /// @param[in]  i   The index of the element to remove
    // ------------------------------------------------------------------------------------------------------
    CUDA_H
    inline void remove_element(const size_t i) 
    {
        // Decrease the number of elements
        --_num_elements;
        
        // First remove the bit which is the element
        _data[i / elements_per_bin].remove_bit(i % elements_per_bin);

        // Now for each of the bins to the right, set the LSB of the container to the MSB of the
        // container to the right, then shift all bits of the container to the left
        for (size_t bin = i / elements_per_bin; bin < _bins; ++bin) {
            _data[bin].set(elements_per_bin - 1, _data[bin + 1].get(0));
            _data[bin + 1].shift_left(1);
        }
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Adds an element to the back of the vector
    /// @param[in]  value     The element to add to the back of the vector
    // ------------------------------------------------------------------------------------------------------
    CUDA_H
    inline void push_back(const byte value)
    {
        // If the container is currently full
        if (!((elements_per_bin * _bins) > _num_elements)) { _data.push_back(internal_container()); ++_bins; }
        
        // Add the new element
        _data[_num_elements / elements_per_bin].set(_num_elements % elements_per_bin, value);
        ++_num_elements;
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts to bool vector
    // ------------------------------------------------------------------------------------------------------
    thrust::host_vector<uint8_t> to_thrust_vector() const 
    {
        standard_container host_vec;
        for (auto i = 0; i < num_elements; ++i)
            host_vec.push_back(get(i));
        return host_vec;
    }
};


}
#endif
