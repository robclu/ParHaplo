// ----------------------------------------------------------------------------------------------------------
/// @file   read_info.hpp
/// @brief  Header file for the ReadInfo class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_READ_INFO_HPP
#define PARAHAPLO_READ_INFO_HPP

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @struct     Read
/// @brief      Holds information about a read, such as the start and end posiitons and the length, it does
///             not hold the read data as this class is intended to be used in conjuncito with a Block. For
///             example when a block has been read (it holds the data) the read can then hold the properties
///             of each of the rows in the Block. The intention is to make the code read better.
// ----------------------------------------------------------------------------------------------------------
class Read {
private:
    int _start;          //!< The start position of the read
    int _end;            //!< The end position of the read
    int _length;         //!< The length of the read 
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Default constructor
    // ------------------------------------------------------------------------------------------------------
    Read() : _start(-1), _end(-1), _length(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Construct a row, setting the start and the end and determining the length
    // ------------------------------------------------------------------------------------------------------
    Read(int start, int end) 
    : _start(start), _end(end), _length(end - start + 1) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Returns the length of the Read
    // ------------------------------------------------------------------------------------------------------
    inline int length() const { return _length; };
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Gets the start index of the read (the first element in a row which is not a gap)
    /// @return THe index of the first non-gap element in a read 
    // ------------------------------------------------------------------------------------------------------
    inline int start() const { return _start; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief  Gets the last index of the read (the last element in the row which is not a gap)
    /// @return The index of the last non-gap element in a read
    // ------------------------------------------------------------------------------------------------------
    inline int end() const { return _end; }
};

}           // End namespace haplo

#endif      // PARAHAPLO_READ_INFO_HPP
