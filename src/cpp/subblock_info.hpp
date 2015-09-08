// ----------------------------------------------------------------------------------------------------------
/// @file   subblock_info.hpp
/// @brief  Header file for the SubBlockInfo class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SUBBLOCK_INFO_HPP
#define PARAHAPLO_SUBBLOCK_INFO_HPP

namespace haplo {
   
// ----------------------------------------------------------------------------------------------------------
/// @struct     SUBBLOCKInfo
/// @brief      Holds information about the unsplittable sub-blocks of another block -- the start column and 
///             the end column of the sub-block within the other block
// ----------------------------------------------------------------------------------------------------------
class SubBlockInfo {
private:
    int _start;         //!< The start index of the sub-block
    int _end;           //!< The end index of the  sub-block
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Default constructor which sets the start and end position to 0
    // ------------------------------------------------------------------------------------------------------
    SubBlockInfo() : _start(0), _end(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Constructs a new SubBlockInfo using the given start and end positions 
    // ------------------------------------------------------------------------------------------------------
    SubBlockInfo(int start, int end) : _start(start), _end(end) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Gets the start index 
    /// @return THe index of the first unsplittable column in the main block
    // ------------------------------------------------------------------------------------------------------
    inline int start() const { return _start; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief  Gets the end index
    /// @return The index of the last unsplittable column in the main block
    // ------------------------------------------------------------------------------------------------------
    inline int end() const { return _end; }
};

}           // End namesapce haplo

#endif      // PARAHAPLO_SUBBLOCK_INFO_HPP
