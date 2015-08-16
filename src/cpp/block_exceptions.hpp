// ----------------------------------------------------------------------------------------------------------
/// @file   block_exceptions.hpp
/// @brief  Header file for all block exceptions for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK_EXCEPTIONS_HPP
#define PARAHAPLO_CPP_BLOCK_EXCEPTIONS_HPP

#include <iostream>
#include <exception>
#include <string>

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class      BlockInputParseError
/// @brief      Defines an exception class for when a Block cannot parse the input file to load the data
// ----------------------------------------------------------------------------------------------------------
class BlockInputParseError : public std::exception {
public:
    std::string _message;                                                   //!< Error message for exception
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Defines the error message for when the exception occurs
    /// @param[in]  filename    The name of the file from which the input coukld not be parsed
    // ------------------------------------------------------------------------------------------------------
    BlockInputParseError(const std::string filename)
    : _message( "Error : Could not parse input file : " +
                filename                                +
                " : correctly to create Block data"     ) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the error message when the exception is caught
    /// @return     The error message for the exception
    // ------------------------------------------------------------------------------------------------------
    const char* what() const throw() { return _message.c_str(); }
};

}       // End namespace haplo

#endif  // PARAHAPLO_CPP_BLOCK_EXCEPTIONS_HPP
