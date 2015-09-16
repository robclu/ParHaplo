// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library -- cpu implementation
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_IMPLEMENTATION_CPU_HPP
#define PARAHAPLO_BLOCK_IMPLEMENTATION_CPU_HPP

#include "block_interface.hpp"

#include <string>
#include <iostream>

namespace haplo {

// Specialization for CPU block implementation    
template <size_t Rows, size_t Cols, size_t Cores>
class BlockImplementation<Rows, Cols, Cores, Device::CPU> {
public:
    // ----------------------------------- TYPEDEFS ---------------------------------------------------------
    using data_container = std::array<haplo::Data, Rows * Cols>;
    // ------------------------------------------------------------------------------------------------------
   
    BlockImplementation() {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of cores available to the block
    /// @return     The number of cores available to the block 
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t num_cores() { return Cores; }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of rows in the block 
    /// @return     The number of rows in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t rows() { return Rows; }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of columns in the block
    /// @return     The number of columns in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t cols() { return Cols; }         
  
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size (numbre of elements) in the block
    /// @return     The number of elements in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t size() { return Rows * Cols; } 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Loads data into the Block from the given file. The function will look for as many elements 
    ///             as are specified by the Block dimensions, so for a Block with 3 rows and 6 columns, the
    ///             function will look for 3 lines each with 6 elements.                                     \n
    ///                                                                                                      \n
    ///             The file is loaded as a memory mapped file since the input data files will likely be     \n
    ///             huge, which should save have a significant performance incease.                          \n
    ///                                                                                                      \n
    ///             See boost memory mapped files for reference:                                             \n
    ///                 http://www.boost.org/doc/libs/1_38_0/libs/iostreams/doc/index.html                   \n
    ///                                                                                                      \n
    ///             Or if you download the boost libraries the source for the memory mapped files is at:     \n
    ///                 boost/iostreams/device/mapped_file.hpp                                               \n
    /// @param[in]  filename        The name of the file to get the input data from.
    // ------------------------------------------------------------------------------------------------------
    void fill(const std::string filename);
};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

template <size_t R, size_t CL, size_t CR>
void BlockImplementation<R, CL, CR, Device::CPU>::fill(const std::string filename)
{
    std::cout << "Getting data ....\n";
}

}           // End namespace haplo

#endif      // PARAHAPLO_BLOCK_CPU_HPP
