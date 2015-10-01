// ----------------------------------------------------------------------------------------------------------
/// @file   variables.hpp
/// @brief  Header file for variables for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_VARIABLES_HPP
#define PARAHAPLO_VARIABLES_HPP

#include <vector>
#include <iostream>

namespace haplo {

class BinaryArray {
public:
    uint8_t _bits;
public: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to initialize the bits to 0
    // ------------------------------------------------------------------------------------------------------
    BinaryArray() : _bits(0) {};
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the bit at position i
    /// @return     The value of the bit at position i
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t operator[](size_t i) const { return (_bits >> i) & 0x01; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the bit at position i
    /// @param[in]  i       THe position of the bit to set
    /// @param[in]  value   The value to set the bit to (0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(size_t i, uint8_t value) 
    { 
        // If the value is valid and values != value, then we must swap the values
        if (value <= 1 && !((_bits >> i & 0x01) && value)) _bits ^= (0x01 << i);
    }
};

// ----------------------------------------------------------------------------------------------------------
/// @class  BnbVariable    
/// @brief  Variable class for the varibles used in the branch and bound implementation for solving the
///         haplotype assembly problem. The variable can be one (linear) or two (non-linear) dimensional
///         variable, and the number of elements in the variables which have been determined is kept track of
// ----------------------------------------------------------------------------------------------------------
template <size_t Dimensions>
class BnbVariable;

// Specialize for the one dimensional case
template <> class BnbVariable<1> {
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor which sets the number of elements
    /// @param[in]  num_elements    The number of elements in the variable
    // ------------------------------------------------------------------------------------------------------
    BnbVariable(const size_t num_elements) 
    : _num_set_elements(0), _data(num_elements / 8 + 1) {}
    
    // Type alias for data container
    using container_type = std::vector<BinaryArray>;
private:
    container_type  _data;              //!< The values of each of the binary variables
    size_t          _num_set_elements;  //!< The number of elements which have been set (determined)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the cariable container
    /// @param[in]  i   The index of the variable in the container to get
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t operator()(const size_t i) const { return _data[i / 8][i % 8]; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const size_t i, const uint8_t value)
    {
        ++_num_set_elements;
        _data[i / 8].set(i % 8, value);       // Not doing error checking here
    }
};

// Specialize for the two dimensional case
template <> class BnbVariable<2> {
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor which sets the number of elements
    /// @param[in]  num_elements_d1     The number of elements in the first dimension 
    /// @param[in]  num_elements_d2     The number of elements in the second dimension 
    // ------------------------------------------------------------------------------------------------------ 
    BnbVariable(const size_t num_elements_d1, const size_t num_elements_d2) 
    : _num_set_elements(0)                              , 
      _num_elements_d1(num_elements_d1)                 ,
      _num_elements_d2(num_elements_d2)                 ,
      _data(num_elements_d1 * num_elements_d2 / 8 + 1)  {}

    // Type alias for container
    using container_type = std::vector<BinaryArray>;
private:
    container_type      _data;              //!< The values of each of the binary variables
    size_t              _num_elements_d1;   //!< The number of elements in the first dimension
    size_t              _num_elements_d2;   //!< The number of elements in the second dimension
    size_t              _num_set_elements;  //!< The number of elements which have been set (determined)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the variable container
    /// @param[in]  i   The index of the variable in the first dimension of the container (row)
    /// @param[in]  j   The index of the variable in the second dimension of the container (column)
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t operator()(const size_t i, const size_t j) const 
    {
        const size_t index = i * _num_elements_d2 + j;
        return _data[index / 8][index % 8];
    }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable in the first dimension of the container (row)
    /// @param[in]  j       The index of the variable in the second dimension of the container (column)
    /// @tparam[in] value   The value to set the variable to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    inline void set(const size_t i, const size_t j, const uint8_t value)
    {
        ++_num_set_elements;
        const size_t index = i * _num_elements_d2 + j;
        _data[index / 8].set(index % 8, value);
    }
};
}           // End namespace haplo

#endif      // PARAHAPLO_VARIABLES_HPP

