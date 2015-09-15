// ----------------------------------------------------------------------------------------------------------
/// @file   branch_and_bound_variables.hpp
/// @brief  Header file for branch and bound variables for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BRANCH_AND_BOUND_VARIABLES_HPP
#define PARAHAPLO_CPP_BRANCH_AND_BOUND_VARIABLES_HPP

#include <vector>
#include <iostream>

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class  BinaryContainer
/// @brief  Class that uses 1 bit per binary variable rather than 8, which, for huge containers, will make a
///         difference -- especially when data must be transferred between the CPU and GPU
// ----------------------------------------------------------------------------------------------------------
class BinaryContainer {
public: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor -- sets all the bits to 0
    // ------------------------------------------------------------------------------------------------------
    BinaryContainer()
    : b0(0), b1(0), b2(0), b3(0), b4(0), b5(0), b6(0), b7(0) {}
        
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of a bit
    /// @param[in]  i       The index of the bit in the container
    /// @param[in]  value   The value to set the bit to 
    // ------------------------------------------------------------------------------------------------------
    void set(uint8_t i, uint8_t value)
    {
        if (value < 2) {
            switch(i) {
                case 0:
                    b0 = value == 0 ? 0 : 1;
                    break;
                case 1:
                    b1 = value == 0 ? 0 : 1;
                    break;
                case 2:
                    b2 = value == 0 ? 0 : 1;
                    break;
                case 3:
                    b3 = value == 0 ? 0 : 1;
                    break;
                case 4:
                    b4 = value == 0 ? 0 : 1;
                    break; 
                case 5:
                    b5 = value == 0 ? 0 : 1;
                    break;
                case 6:
                    b6 = value == 0 ? 0 : 1;
                    break; 
                case 7:
                    b7 = value == 0 ? 0 : 1;
                    break; 
                default: break;
            }
        }
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Get a value of a bit
    /// @param[in]  i   The index of the bit to get
    /// @return     The value of the bit i 
    // ------------------------------------------------------------------------------------------------------
    uint8_t get(uint8_t i) const
    {
        switch(i) {
            case 0:
                return b0;
            case 1:
                return b1;
            case 2:
                return b2;
            case 3:
                return b3; 
            case 4:
                return b4;
            case 5:
                return b5; 
            case 6:
                return b6;
            case 7:
                return b7; 
            default: break;
        }
    }   
private:
    uint8_t b0 : 1;         //!< Zeroth bit
    uint8_t b1 : 1;         //!< First bit
    uint8_t b2 : 1;         //!< Second bit
    uint8_t b3 : 1;         //!< Third bit
    uint8_t b4 : 1;         //!< Forth bit
    uint8_t b5 : 1;         //!< Fifth bit
    uint8_t b6 : 1;         //!< Sixth bit
    uint8_t b7 : 1;         //!< Seventh bit 
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
    using container_type = std::vector<BinaryContainer>;
private:
    container_type  _data;              //!< The values of each of the binary variables
    size_t          _num_set_elements;  //!< The number of elements which have been set (determined)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a value from the cariable container
    /// @param[in]  i   The index of the variable in the container to get
    // ------------------------------------------------------------------------------------------------------
    uint8_t get(const size_t i) const { return _data[i / 8].get(i % 8); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable to set
    /// @param[in]  value   The value to set the varibale to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    void set(const size_t i, const uint8_t value)
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
    using container_type = std::vector<BinaryContainer>;
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
    uint8_t get(const size_t i, const size_t j) const 
    {
        const size_t index = i * _num_elements_d2 + j;
        return _data[index / 8].get(index % 8);
    }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a value in the variable container
    /// @param[in]  i       The index of the variable in the first dimension of the container (row)
    /// @param[in]  j       The index of the variable in the second dimension of the container (column)
    /// @tparam[in] value   The value to set the variable to (must be 0 or 1)
    // ------------------------------------------------------------------------------------------------------
    void set(const size_t i, const size_t j, const uint8_t value)
    {
        ++_num_set_elements;
        const size_t index = i * _num_elements_d2 + j;
        _data[index / 8].set( index % 8, value);
    }
};
}           // End namespace haplo

#endif      // PARAHAPLO_CPP_BRANCH_AND_BOUND_VARIABLES_HPP

