// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo link class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_LINK_HPP
#define PARHAPLO_LINK_HPP

#include "devices.hpp"

#include <tbb/tbb.h>
#include <memory>

namespace haplo {

static constexpr uint8_t small = 0x00;
static constexpr uint8_t big   = 0x01;

// ----------------------------------------------------------------------------------------------------------
/// @class      Link
/// @brief      A link between two nodes, there is a homozygous component -- how strongly correlated the nodes
///             are (that they should have the same value) -- and a heterozygous component -- how stongly they
///             should be different.
// ----------------------------------------------------------------------------------------------------------
class Link {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using atomic_type = tbb::atomic<size_t>;
private:  
    atomic_type     _homo_weight;        //!< Weight of the link if the nodes have the same ideal values
    atomic_type     _hetro_weight;       //!< Weight of the link if the nodes have different ideal values
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Link() noexcept : _homo_weight{0}, _hetro_weight{0} {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor with arguments
    // ------------------------------------------------------------------------------------------------------
    Link(atomic_type homo_weight, atomic_type hetro_weight) noexcept 
    : _homo_weight(homo_weight), _hetro_weight(hetro_weight) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for link class
    // ------------------------------------------------------------------------------------------------------
    ~Link() noexcept {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the homozygous weight 
    /// @return     A reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& homo_weight() { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the heteroygous weight 
    /// @return     A reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& hetro_weight() { return _hetro_weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the homozygous weight 
    /// @return     A tosnt reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& homo_weight() const { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the heteroygous weight 
    /// @return     A const reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& hetro_weight() const { return _hetro_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the link so that it cant be used bya sorting function
    /// @return     The value (maximum weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return std::max(_homo_weight, _hetro_weight); }
};


// ----------------------------------------------------------------------------------------------------------
/// @class      LinkContainer 
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class LinkContainer;

}           // End namespace haplo
#endif      // PARAHAPLO_LINK_HPP
