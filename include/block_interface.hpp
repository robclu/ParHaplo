// ----------------------------------------------------------------------------------------------------------
/// @file   block_interface.hpp
/// @brief  Header file for the block interface class for the parahaplo library -- this just defines an 'empty' 
///         base which allows the specific block implementation to be selected at compile time 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_INTERFACE_HPP
#define PARAHAPLO_BLOCK_INTERFACE_HPP

#include "data.hpp"
#include "devices.hpp"

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class  BlockInterface  
/// @brief  Empty base class for Blocks which provides an interface for selecting the dvice type to use for
///         all block operations -- CPU, GPU ot PHI
/// @tparam Implementation  The type of implementation to use -- CPU, GPU, PHI
// ----------------------------------------------------------------------------------------------------------
template <typename Implementation>
class BlockInterface : Implementation {
public:
    // Type alias for the block type 
    using block_type = Implementation;

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a pointer to the base class (block implementation)
    /// @return     A pointer to the base class
    // ------------------------------------------------------------------------------------------------------
    Implementation* implementation() { return static_cast<Implementation*>(*this); }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a contant pointer to the base class (implementation)
    /// @return     A constant pointer to the base class
    // ------------------------------------------------------------------------------------------------------
    const Implementation* implementation() const { return static_cast<const Implementation*>(this); } 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of cores that are available to the block
    /// @return     The number of cores available to the block
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t num_cores() { return block_type::num_cores(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the block 
    /// @return     The size of the block 
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t size() { return block_type::size(); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data
    /// @param[in]  filename    The name of the file which containes the data to fill the block with
    // ------------------------------------------------------------------------------------------------------
    void fill(const std::string filename) { implementation()->fill(filename); }
};

// ----------------------------------------------------------------------------------------------------------
/// @class  BlockImplementation
/// @brief  Block implementation class which can be specialized to provide CPU, GPU and PHI implementations 
/// @tparam Rows        The number of rows in the block
/// @tparam Cols        The number of columns in the block
/// @tparam Cores       The number of cores available to the block -- later we can change this to be a custom 
///         data type which allows cores for each dimension to be specified
/// @tparam DeviceType  The type of device the block uses -- CPU, GPU, or PHI 
// ----------------------------------------------------------------------------------------------------------
template <size_t Rows, size_t Cols, size_t Cores, haplo::Device DeviceType>
class BlockImplementation;

}               // End namespace haplo

#endif          // PARAHAPLO_BLOCK_INTERFACE_HPP
