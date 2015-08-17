// ----------------------------------------------------------------------------------------------------------
/// @file   data.hpp
/// @brief  Header file for the parahaplo data class which defines input data
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_DATA_HPP
#define PARAHAPLO_DATA_HPP

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class  Data 
/// @brief  Holds a data input for the haplotype assembly problem. Since inputs can only be an element of {0,
///         1, -}, we only need 2 bits to represent each data, which will save space 
// ----------------------------------------------------------------------------------------------------------
class Data {
private:
    unsigned int _value : 2;                                                    //!< Value of data {0, 1, -}

public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Takes a character as input (sicne that's what the input will be from the input file, and
    ///             converts it to the relevant 2 bit value
    /// @param[in]  value   The char value of the data
    // ------------------------------------------------------------------------------------------------------
    Data(const char value);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the data
    /// @return     The value of the data
    // ------------------------------------------------------------------------------------------------------
    inline unsigned int value() const { return _value; }
};

// ------------------------------------------ Implementation ------------------------------------------------

Data::Data(const char value) 
{
    switch(value) {
        case '0':
            _value = 0;
            break;
        case '1':
            _value = 1;
            break;
        case '-':
            _value = 2;
            break;
        default:            // Anything else is invalid
            _value = 3;
    }    
}

}       // End namespace haplo

#endif // PARAHAPLO_DATA_HPP
