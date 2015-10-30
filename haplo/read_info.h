// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo read info class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_READ_INFO_H
#define PARHAPLO_READ_INFO_H

#include "cuda_defs.h"

namespace haplo {
           
// ----------------------------------------------------------------------------------------------------------
/// @class      ReadInfo
/// @brief      Stores some information about the reads 
// ----------------------------------------------------------------------------------------------------------
class ReadInfo {
private:
    size_t      _start_idx;
    size_t      _end_idx;
    size_t      _offset;
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    ReadInfo(){}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the values
    /// @param[in]  read_number     The number of the read (row number)
    /// @param[in]  start_index     The start index (column) of the read
    /// @param[in]  end_index       The end index of the read
    /// @param[in]  offset          The offset in memory of the start of the read
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    ReadInfo(const size_t a, const size_t start_idx, const size_t end_idx, const size_t offset   )
    : _start_idx(start_idx), _end_idx(end_idx), _offset(offset) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the staet index of the read 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t start_index() const { return _start_idx; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the start index of the read
    /// @param[in]  value   The new value for the start index
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void set_start_index(const size_t value) { _start_idx = value; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the end index of the read
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t end_index() const { return _end_idx; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the end index of the read
    /// @param[in]  value   The new value for the end index
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void set_end_index(const size_t value) { _end_idx = value; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the offset of the read in memoy
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t offset() const { return _offset; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the length of the read 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t length() const { return _end_idx - _start_idx + 1; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns true if the element is in the read
    /// @param[in]  index   The index of the element
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline bool element_exists(const size_t index) const { return (index >= _start_idx && index <= _end_idx); }
};

}           // End namespace haplo
#endif      // PARAHAPLO_READ_INFO_H
