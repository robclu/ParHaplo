// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo link class for use with GPUs
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_LINK_GPU_CUH
#define PARHAPLO_LINK_GPU_CUH

#ifdef __CUDACC__
#define CUDA_HD __host__ __device__
#else
#define CUDA_HD
#endif

#include <cuda.h>
#include <math.h>

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class      Link
/// @brief      A link between two nodes, there is a homozygous component -- how strongly correlated the nodes
///             are (that they should have the same value) -- and a heterozygous component -- how stongly they
///             should be different.
// ----------------------------------------------------------------------------------------------------------
class LinkGpu {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
private:  
    size_t     _homo_weight;        //!< Weight of the link if the nodes have the same ideal values
    size_t     _hetro_weight;       //!< Weight of the link if the nodes have different ideal values
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    LinkGpu() noexcept : _homo_weight{0}, _hetro_weight{0} {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor with arguments
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    LinkGpu(const size_t homo_weight, const size_t hetro_weight) noexcept 
    : _homo_weight(homo_weight), _hetro_weight(hetro_weight) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for link class
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    ~LinkGpu() noexcept {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the homozygous weight 
    /// @return     A reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t& homo_weight() { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the heteroygous weight 
    /// @return     A reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t& hetro_weight() { return _hetro_weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the homozygous weight 
    /// @return     A tosnt reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline const size_t& homo_weight() const { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the heteroygous weight 
    /// @return     A const reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline const size_t& hetro_weight() const { return _hetro_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the link so that it cant be used bya sorting function
    /// @return     The value (maximum weight of the node)
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t max() const 
    {
        #ifdef __CUDAACC__
            return max(_homo_weight, _hetro_weight); 
        #else
            return std::max(_homo_weight, _hetro_weight); 
        #endif         
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the min value of the link
    /// @return     The min value of the link
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t min() const 
    {
        #ifdef __CUDAACC__
            return min(_homo_weight, _hetro_weight); 
        #else
            return std::min(_homo_weight, _hetro_weight); 
        #endif 
    }
};


// ----------------------------------------------------------------------------------------------------------
/// @class      LinkContainer 
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class LinkContainer;

}           // End namespace haplo
#endif      // PARAHAPLO_LINK_HPP
