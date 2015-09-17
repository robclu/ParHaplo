// ----------------------------------------------------------------------------------------------------------
/// @file   usplit_block_cpu.hpp
/// @brief  Header file for the unsplittable block class cpu implementation for the parahaplo library 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_USPLIT_BLOCK_IMPLEMENTATION_CPU_HPP
#define PARAHAPLO_USPLIT_BLOCK_IMPLEMENTATION_CPU_HPP

#include "block_interface.hpp"
#include "usplit_block_interface.hpp"

#include <tbb/concurrent_unordered_map.h>
#include <bitset>
#include <numeric>

namespace haplo {

// Specialization for using the CPU
template <typename BaseBlock>
class UnsplittableBlockImplementation<Device::CPU, BaseBlock> : BaseBlock {
public:
    // --------------------------------------- ALIAS'S ------------------------------------------------------
    using base_type         = BaseBlock;
    using data_container    = std::vector<typename base_type::data_type>;
    using subinfo_type      = typename base_type::subinfo_type;
    using ston_container    = typename std::array<bool, base_type::rows()>;
private:
    data_container  _data;              //!< Data for the unsplittable block
    ston_container  _singleton_info;    //!< Which rows are singletons and which are not
    size_t          _index;             //!< The index of the unsplittable block in the base block
    size_t          _size;              //!< The size -- total number of elements -- in the unsplittable block
    size_t          _rows;              //!< The number of rows in the unsplittable block
    size_t          _cols;              //!< The number of columns in the unsplitatble block
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor which sets the base block for the unsplittable block 
    /// @param[in]  base_block  The block which is being used as the base of this class. This unsplittable
    ///             block is a sub-region of the base block, 
    /// @param[in]  index       The index of this unsplittable block in the base block
    // ------------------------------------------------------------------------------------------------------
    UnsplittableBlockImplementation(base_type& base_block, const size_t index = 0);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the base block 
    /// @return     A pointer to the base block 
    // ------------------------------------------------------------------------------------------------------
    base_type* base_block() { return static_cast<base_type*>(this); }
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      First removes all the singleton rows (rows with only 1 value) from the unsplittable block
    ///             from which the number of rows, and hence the total size of the block, can be determined.
    // ------------------------------------------------------------------------------------------------------
    void determine_params(); 
};

// -------------------------------------- IMPLEMENTATIONS ---------------------------------------------------

// ------------------------------------------ PUBLIC --------------------------------------------------------

template <typename BaseBlock>
UnsplittableBlockImplementation<Device::CPU, BaseBlock>::UnsplittableBlockImplementation(base_type& base_block, 
                                                                                         const size_t index)
: base_type(base_block), _data(0), _index(index) 
{
    determine_params();
    std::cout << " Working!";
}

template <typename BaseBlock>
void UnsplittableBlockImplementation<Device::CPU, BaseBlock>::determine_params()
{
    // Number of rows in the Expression 
    constexpr size_t rows           = base_type::rows();
    constexpr size_t num_threads    = base_type::num_cores();
    
    // Create threads where each checks for singletons
    const size_t threads_to_use = num_threads > rows ? rows : num_threads;

    // Information for the sub-block
    const subinfo_type& info = base_block()->subblock_info(_index);
        
    // Check which rows are singletons
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices) 
        {
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                size_t thread_iters = ops::get_thread_iterations(idx, rows, threads_to_use);
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t non_gaps     = 0;                        // Number of row elements which aren't gaps
                    size_t row_offset   = ops::thread_map(idx, threads_to_use, it);
                    // Now we need to go through all elements in the row 
                    // of base block's rows which are part of this unspittable 
                    // block and check of there is only a single element
                    for (size_t col_offset = info.start(); col_offset <= info.end() && non_gaps < 2; ++col_offset) {
                        if (base_block()->operator()(row_offset, col_offset).value() != 2) 
                            ++non_gaps;                             // Not a gap, so increment num_elements
                    }
                    // If we have found more than a single element then the row needs to be added
                    // to the UnsplittableBlock, so set the relevant value to 1 (not a singleton)
                    _singleton_info[row_offset] = non_gaps > 1 ? 1 : 0;
                }
            }
        }
    );
    
    // This acclally counts the number of rows which aren't singletons
    _rows = std::accumulate(_singleton_info.begin(), _singleton_info.end(), 0);
    _cols = info.columns();
    _size = _rows * info.columns();
    _data.resize(_size, static_cast<uint8_t>(0));   
}

}               // End namespace haplo

#endif          // PARAHAPLO_USPLIT_BLOCK_IMPLEMENTATION_CPU_HPP
