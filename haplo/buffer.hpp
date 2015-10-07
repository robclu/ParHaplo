// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo buffer container 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_BUFFER_HPP
#define PARHAPLO_BUFFER_HPP

namespace haplo {
           
// ----------------------------------------------------------------------------------------------------------
/// @class      Buffer
/// @brief      Wrapper for a buffer object
/// @tparam     Type The type of object to create a buffer of
// ----------------------------------------------------------------------------------------------------------
template <typename Type>
class Buffer {
private:
    void* _ptr;         //!< The pointer for the buffer
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- tries to create a buffer with a given size
    /// @param[in]  num_bytes   The number of bytes for the buffer
    // ------------------------------------------------------------------------------------------------------
    Buffer(const size_t num_bytes) : _ptr(operator new(num_bytes, std::nothrow)) {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for the buffer
    // ------------------------------------------------------------------------------------------------------
    ~Buffer() { operator delete(_ptr); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      For determining of the buffer was created successfully
    /// @return     If the buffer was created successfully or not
    // ------------------------------------------------------------------------------------------------------
    operator bool() const { return _ptr; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fetches the buffer data
    /// @return     The pointer to the buffer data
    // ------------------------------------------------------------------------------------------------------
    Type* fetch_data() const { return static_cast<Type*>(_ptr); }
};

}           // End namespace haplo
#endif      // PARAHAPLO_BUFFER_HPP
