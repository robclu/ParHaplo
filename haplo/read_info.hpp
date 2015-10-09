// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo read info class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_READ_INFO_HPP
#define PARHAPLO_READ_INFO_HPP

namespace haplo {
           
// ----------------------------------------------------------------------------------------------------------
/// @class      ReadInfo
/// @brief      Stores some informatio about each of the reads
// ----------------------------------------------------------------------------------------------------------
class ReadInfo {
private:
    size_t      _read_number;       
    size_t      _start_idx;
    size_t      _end_idx;
    size_t      _offset;
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    ReadInfo() noexcept {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the values
    /// @param[in]  read_number     The number of the read (row number)
    /// @param[in]  start_index     The start index (column) of the read
    /// @param[in]  end_index       The end index of the read
    /// @param[in]  offset          The offset in memory of the start of the read
    // ------------------------------------------------------------------------------------------------------
    ReadInfo(const size_t read_number, const size_t start_idx, 
             const size_t end_idx    , const size_t offset   ) noexcept
    : _read_number(read_number), _start_idx(start_idx), _end_idx(end_idx), _offset(offset) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the read number of the read information
    // ------------------------------------------------------------------------------------------------------
    inline size_t read_number() const { return _read_number; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the staet index of the read 
    // ------------------------------------------------------------------------------------------------------
    inline size_t start_index() const { return _start_idx; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the start index of the read
    /// @param[in]  value   The new value for the start index
    // ------------------------------------------------------------------------------------------------------
    inline void set_start_index(const size_t value) { _start_idx = value; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the end index of the read
    // ------------------------------------------------------------------------------------------------------
    inline size_t end_index() const { return _end_idx; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the end index of the read
    /// @param[in]  value   The new value for the end index
    // ------------------------------------------------------------------------------------------------------
    inline void set_end_index(const size_t value) { _end_idx = value; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the offset of the read in memoy
    // ------------------------------------------------------------------------------------------------------
    inline size_t offset() const { return _offset; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the length of the read 
    // ------------------------------------------------------------------------------------------------------
    inline size_t length() const { return _end_idx - _start_idx + 1; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns true if the element is in the read
    /// @param[in]  index   The index of the element
    // ------------------------------------------------------------------------------------------------------
    inline bool element_exists(const size_t index) const { return (index >= _start_idx && index <= _end_idx); }
};

}           // End namespace haplo
#endif      // PARAHAPLO_READ_INFO_HPP
