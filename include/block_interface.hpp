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
    // ------------------------------------------ ALIASES ---------------------------------------------------
    using block_type        = Implementation;
    using reference_type    = typename block_type::reference_type;
    using data_type         = typename block_type::data_type;
    using read_type         = typename block_type::read_type;
    using subinfo_type      = typename block_type::subinfo_type;
    // ------------------------------------------------------------------------------------------------------
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for the block, this constructor must be used otherwise the behaviour is
    ///             undefined since the block will not have any data
    /// @param[in]  data_file   The name of the data file to get the data from
    // ------------------------------------------------------------------------------------------------------
    explicit BlockInterface(const char* data_file) : Implementation(data_file) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a pointer to the base class (block implementation)
    /// @return     A pointer to the base class
    // ------------------------------------------------------------------------------------------------------
    Implementation* implementation() { return static_cast<Implementation*>(this); }
   
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

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of rows in the block 
    /// @return     The number of rows in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t rows() { return block_type::rows(); }
   
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of columns in the block
    /// @return     The number of columns in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t cols() { return block_type::cols(); }     

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the block 
    /// @return     The size of the block 
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t size() { return block_type::size(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the information for the sub block i of this block. This will throw an out of range
    ///             error if the sub-block for which the information is wanted is not valid
    /// @param[in]  i   The index of the sub-block for which the information must be given
    /// @return     The subblock information for the sub-block with index i
    // ------------------------------------------------------------------------------------------------------
    const subinfo_type& subblock_info(const size_t i) const { return implementation()->subblock_info(i); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a refernce to the element of the block at row and column -- does not implement bound
    ///             checking as this is optimozed for performance, .at() can be implemented at a later stage 
    ///             if bound checking access is required
    /// @param[in]  row     The row of the element in the block
    /// @param[in]  col     The column of the element in the block 
    // ------------------------------------------------------------------------------------------------------
    reference_type operator()(size_t row, size_t col) { return implementation()->operator()(row, col); }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the read at row i of the block -- nor error checking for performance, again .at can
    ///             be implemented at a later stage
    /// @param[in]  The index of the read to get from the block
    /// @return     The read at row i in the block
    // ------------------------------------------------------------------------------------------------------
    const read_type& operator[](size_t i) const { return implementation()->operator[](i); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data
    /// @param[in]  filename    The name of the file which containes the data to fill the block with
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* filename) { implementation()->fill(filename); }
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
