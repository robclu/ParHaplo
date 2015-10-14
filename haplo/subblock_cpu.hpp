/// @file   subblock_cpu.hpp
/// @brief  Header file for the subblock class cpu implementation for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SUB_BLOCK_CPU_HPP
#define PARAHAPLO_SUB_BLOCK_CPU_HPP

#include "devices.hpp"
#include "processor_cpu.hpp"
#include "tree_cpu.hpp"
#include "subblock.hpp"
#include "snp_info_gpu.h"

#include <sstream>

namespace haplo {

// Specialization for the CPU implementation of the unsplittable block
template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY>
class SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu> : public BaseBlock {
public:
    // ------------------------------------------- ALIAS'S --------------------------------------------------`
    using sub_block_type        = SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>;
    using atomic_type           = tbb::atomic<size_t>;
    using tree_type             = Tree<sub_block_type, devices::cpu>;
    using binary_vector         = BinaryVector<2>;              
    using atomic_vector         = tbb::concurrent_vector<size_t>;
    using node_container        = NodeContainer<devices::cpu>;
    using concurrent_umap       = typename BaseBlock::concurrent_umap;
    using read_info_container   = typename BaseBlock::read_info_container;
    using snp_info_container    = typename BaseBlock::snp_info_container;
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t     THREADS_X   = ThreadsX;
    static constexpr size_t     THREADS_Y   = ThreadsY;
private:
    size_t              _num_nih;           //!< The number of NIH columns
    size_t              _index;             //!< The index of the unsplittable block within the base block
    size_t              _cols;              //!< The number of columns in the sub block
    size_t              _rows;              //!< The number of rows in the sub block 
    size_t              _elements;          //!< The number of elements in the sub block
    size_t              _base_start_row;    //!< The start row of the subblock in the base block
    
    binary_vector       _data;              //!< The data for the block
    binary_vector       _haplo_one;         //!< The first haplotype
    binary_vector       _haplo_two;         //!< The second haplotype 
    binary_vector       _alignments;        //!< The alignments of the haplotypes         
    tree_type           _tree;              //!< The tree to be solved for the block
        
    read_info_container _read_info;         //!< The information for each of the reads (rows)
    snp_info_container  _snp_info;          //!< The information for each of the snps (columns)

    // These variables are for making the processing faster
    concurrent_umap     _duplicate_rows;        //!< Map of duplicate rows 
    concurrent_umap     _duplicate_cols;        //!< Map of duplicate cols
    concurrent_umap     _row_multiplicities;    //!< How many duplicates each row has
    
    // Friend class that can process rows and columns    
    template <typename FriendType, byte ProcessType, byte DeviceType>
    friend class Processor;
    
    // Tree is s friend class so that it can access the data 
    template <typename SubBlockType, byte DeviceType>
    friend class Tree;
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for when the size (number of elements) is not given (this is the preferred way
    ///             as the _data can be built minimally -- i,e we go through the data from the base block
    ///             which makes up this block and then add only non-singluar rows). Note: Resizing the data
    ///             container (to hold more elements) is not expensive -- adding 16 elements is the equivalent 
    ///             of creating a single int
    /// @param[in]  block   The block from which this block derives
    /// @param[in]  index   The index of the unsplittable block within blokc (block has a specific number of 
    ///             unsplittable blocks which can be made from it)
    // ------------------------------------------------------------------------------------------------------
    explicit SubBlock(const BaseBlock& block, const size_t index);
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the element at position row_idx, col_idx
    /// @param[in]  row_idx     The index of the row of the element
    /// @param[in]  col_idx     The index of the column of the element
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t operator()(const size_t row_idx, const size_t col_idx) const ;
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the sub block
    // ------------------------------------------------------------------------------------------------------
    inline size_t index() const { return _index; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of reads that make up the sub block
    // ------------------------------------------------------------------------------------------------------
    inline size_t reads() const { return _rows; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a constant reference to the tree for the block
    /// @return     A constant ference to the tree for the block
    // ------------------------------------------------------------------------------------------------------
    inline const tree_type& tree() const { return _tree; }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Searches the tree to find the haplotypes
    // ------------------------------------------------------------------------------------------------------
    inline void find_haplotypes() { _tree.template explore<ThreadsX, ThreadsY>(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the number of elements in the subblock
    // ------------------------------------------------------------------------------------------------------
    inline size_t size() const { return _elements; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Prints out the haplotypes for the sub-block
    // ------------------------------------------------------------------------------------------------------
    void print_haplotypes() const;

    // ------------------------------------------------------------------------------------------------------
    /// @brief      A refernce to the data
    // ------------------------------------------------------------------------------------------------------
    inline binary_vector& data()  { return _data; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      A reference to the first haplotype
    // ------------------------------------------------------------------------------------------------------
    inline const binary_vector& haplo_one() const { return _haplo_one; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      A reference to the second haplotype 
    // ------------------------------------------------------------------------------------------------------
    inline const binary_vector& haplo_two() const { return _haplo_two; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      A reference to the alignments
    // ------------------------------------------------------------------------------------------------------
    inline const binary_vector& alignments() const { return _alignments; }
 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the start row of the sub block in the base block 
    // ------------------------------------------------------------------------------------------------------
    inline size_t base_start_row() const { return _rows; }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the MEC score of the haplotpye 
    // ------------------------------------------------------------------------------------------------------
    void determine_mec_score() const;

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a reference to the read information
    // ------------------------------------------------------------------------------------------------------
    inline thrust::host_vector<ReadInfo>& read_info() { return _read_info; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the snp info as a host vector
    // ------------------------------------------------------------------------------------------------------
    inline thrust::host_vector<SnpInfoGpu> snp_info() const 
    {
        thrust::host_vector<SnpInfoGpu> host_snps;
        // Move snps from hash table to vector
        for (auto i = 0; i < _cols; ++i) {
            if (_snp_info.find(i) != _snp_info.end())
                host_snps.push_back(_snp_info.at(i));
            else 
                std::cerr << "error!\n";
        }
        return host_snps;
    }

    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------


    // ------------------------------------------------------------------------------------------------------
    /// @brief      prints the subblock
    // ------------------------------------------------------------------------------------------------------
    void print() const 
    {
        for (size_t r = 0; r < _rows; ++r) {
            for (size_t c = 0; c < _cols; ++c) 
                std::cout << static_cast<unsigned>(operator()(r, c)) << " ";
            std::cout << "\n";
        }
    }

private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a pointer to the block which is a base class of this unsplittable block
    /// @return     A pointer the the block which is a base class of this unsplittable block
    // ------------------------------------------------------------------------------------------------------
    const BaseBlock* base_block() const { return static_cast<const BaseBlock*>(this); }

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the start column index of the subblock in the base block
    // ------------------------------------------------------------------------------------------------------
    inline size_t base_start_index() const { return base_block()->subblock(_index); }

    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the end column index of the subblock in the base block
    //-------------------------------------------------------------------------------------------------------
    inline size_t base_end_index() const { return base_block()->subblock(_index + 1); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the links between the nodes representing the haplotype solution
    // ------------------------------------------------------------------------------------------------------
    void determine_haplo_links();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the data for the unsplittable block with the releavant data from the base block
    // ------------------------------------------------------------------------------------------------------
    void fill();

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Find the duplicate rows
    // ------------------------------------------------------------------------------------------------------
    void find_duplicate_rows();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the parameters for a column -- the start and end index
    /// @param[in]  col_idx     The index of the column to set the parameters for
    /// @param[in]  row_idx     The index of the row to update in teh column info
    /// @param[in]  value       The value of the element at row_idx, col_idx
    // ------------------------------------------------------------------------------------------------------
    void set_col_params(const size_t col_idx, const size_t row_idx, const uint8_t value);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Adds elements to the data vector and updates the offset
    /// @param[in]  row_idx         The index of the row that the elements are being added for
    /// @param[in]  read_length     The number of elements to add from the read
    /// @param[in]  mono_weights    The weights for how many monotone columes are before the start index
    /// @param[in]  offset          The offset in memory for where to start adding the elements
    // ------------------------------------------------------------------------------------------------------
    size_t add_elements(const size_t               row_idx     , const size_t read_length, 
                        const std::vector<size_t>& mono_weights, size_t       offset     );    
};

// -------------------------------------------- IMPLEMENTATIONS ---------------------------------------------

// ----------------------------------------------- PUBLIC ---------------------------------------------------

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY>
SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::SubBlock(const BaseBlock& block, 
                                                                const size_t     index) 
: BaseBlock(block)                                                      , 
  _num_nih(0)                                                           ,
  _index(index)                                                         , 
  _cols(block.subblock(index + 1) - block.subblock(index) + 1)          ,
  _rows(0)                                                              ,                        
  _elements(0)                                                          ,
  _base_start_row(0)                                                    ,
  _data(0)                                                              ,
  _tree(*this, block.subblock(index + 1) - block.subblock(index) + 1)   ,
  _read_info(0)                                                 
{
    std::ostringstream error_message;
    error_message   << "Index for unsplittable block past max index\n" 
                    << "\tindex : "     << index 
                    << "\tmax_index : " << block.num_subblocks() - 1 << "\n"; 
    
    // Check that the index is valid 
    try {
        if (index > block.num_subblocks() - 1) 
             throw std::out_of_range(error_message.str());
    } catch (const std::out_of_range& oor) { 
         std::cerr << "Out of Range error: " << oor.what() << '\n'; 
    }
    
    fill();                                             // Fill the block with data
    _tree.resize(_cols);                                // Resize the tree incase there were monotones
    find_duplicate_rows();                              // Find the duplicate rows and the row mltiplicities
    determine_haplo_links();                            // Find the links between haplotype positions
    _haplo_one.resize(_cols);                           // Allocate memory for haplo one
    _haplo_two.resize(_cols);                           // Allocate memory for haplo two
    _alignments.resize(_rows);                          // Resize the aligments vector
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY> 
uint8_t SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::operator()(const size_t row_idx,
                                                                          const size_t col_idx) const
{
    // If the element exists
    return _read_info[row_idx].element_exists(col_idx) == true 
        ? _data.get(_read_info[row_idx].offset() + col_idx - _read_info[row_idx].start_index()) : 0x03; 
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY> 
void SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::print_haplotypes() const 
{
    for (auto i = 0; i < _haplo_one.size() + 6; ++i) std::cout << "-";
    std::cout << "\nh  : "; 
    for (auto i = 0; i < _haplo_one.size(); ++i) std::cout << static_cast<unsigned>(_haplo_one.get(i));
    std::cout << "\nh` : ";
    for (auto i = 0; i < _haplo_two.size(); ++i) std::cout << static_cast<unsigned>(_haplo_two.get(i));
    std::cout << "\n";
    for (auto i = 0; i < _haplo_two.size() + 6; ++i) std::cout << "-";
}


// -------------------------------------------- PRIVATE -----------------------------------------------------

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY> 
void SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::fill()
{
    size_t offset = 0; size_t monos_found = 0; bool first_row_set = false;
    std::vector<size_t> mono_weights(base_end_index() - base_start_index() + 1);
    
    for (size_t col_idx = base_start_index(); col_idx <= base_end_index(); ++col_idx) {
        if (base_block()->is_monotone(col_idx)) ++monos_found;
        mono_weights[col_idx - base_start_index()] = monos_found;
    }
    _cols -= monos_found;       // Subtract the number of montone columns from the total columns
    
    // Go over each of the data rows and check for singularity
    for (size_t row_idx = 0; row_idx < base_block()->reads(); ++row_idx) {
        // If the row is part of this subblock
        if (base_block()->read_info(row_idx).start_index() >= base_start_index() &&
            base_block()->read_info(row_idx).end_index()   <= base_end_index()    ) {

            // Determine the parameters of the read
            auto read_length = base_block()->read_info(row_idx).length();

            // If the read is not singular
            if (read_length > 1) {
                _read_info.push_back(ReadInfo(_rows, 0, 0, offset));
                offset = add_elements(row_idx, read_length, mono_weights, offset);
                _elements += _read_info[_rows].length();
                ++_rows;
                
                // Check if we found the first row
                if (!first_row_set && offset > 0) 
                    first_row_set = true;
                else if (!first_row_set)
                    ++_base_start_row;
            }
        }
    }
   
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY>
size_t SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::add_elements(
                                                            const size_t               base_row_idx,
                                                            const size_t               read_length ,
                                                            const std::vector<size_t>& mono_weights,
                                                            size_t                     offset      )
{
    // Make sure there is enough space
    _data.resize(_data.size() + read_length);                       
 
    // The start and end index of the column
    size_t read_start     = base_block()->read_info(base_row_idx).start_index() - base_start_index();
    size_t start_col      = read_start - mono_weights[read_start];
    size_t end_col        = start_col + read_length;
    
    // For the case where the start element of a read is monotone, re add the mono weight
    if (base_block()->is_monotone(read_start + base_start_index())) {
        start_col = read_start; end_col = start_col + read_length;
    }
    
    bool   start_set    = false;        // If the start element has been found
    size_t num_elements = 0;            // Number of elements in the read
    size_t mono_counter = 0;            // Number of monotones as the start of the read
    
    for (size_t col_idx = start_col; col_idx < end_col; ++col_idx) {
        auto   base_col_idx  = col_idx + base_start_index() + mono_weights[read_start];
        auto   base_elem_val = base_block()->operator()(base_row_idx, base_col_idx);
        bool   is_mono_col   = base_block()->is_monotone(base_col_idx);

        // If not a monotone column, and start is not found, set start
        // otherwise count the number of initial monotone columns
        if (!start_set && !is_mono_col) { 
            _read_info[_rows].set_start_index(col_idx - mono_counter);
            start_set = true;
        } else if (!start_set && is_mono_col) ++mono_counter;
        
        // Check to see if the column is NIH
        if (!is_mono_col && !base_block()->is_intrin_hetro(base_col_idx)) 
            _snp_info[col_idx].set_type(NIH);
    
        // Check what value to add to the data
        if (base_elem_val == 0 && !is_mono_col) {
            _data.set(offset++, ZERO);
            set_col_params(col_idx, _rows, ZERO);
            ++num_elements;
        } else if (base_elem_val == 1 && !is_mono_col) {
            _data.set(offset++, ONE);
            set_col_params(col_idx, _rows, ONE);
            ++num_elements;
        } else if (base_elem_val == 2 && !is_mono_col) {
            _data.set(offset++, TWO);
            ++num_elements;
        }
    }
    // Set the end index
    _read_info[_rows].set_end_index(_read_info[_rows].start_index() + num_elements - 1);
    
    return offset;
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY>
void SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::set_col_params(const size_t   col_idx,
                                                                           const size_t   row_idx,
                                                                           const uint8_t  value  )
{
    if (_snp_info.find(col_idx) == _snp_info.end()) {
        // Not in map, so set start index to row index
        _snp_info[col_idx] = SnpInfo(row_idx, row_idx);
    } else {
        // In map, so start is set, set end 
        _snp_info[col_idx].end_index() = row_idx;
    }
    // Update the value counter
    value == ZERO ? _snp_info[col_idx].zeros()++
                  : _snp_info[col_idx].ones()++;
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY>
void SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::find_duplicate_rows()
{
    // Create a processor for the rows to determine duplicates
    Processor<sub_block_type, proc::row_dups, devices::cpu> row_processor(*this);
    
    // For each of the rows, from back to front
    for (size_t row_idx = _rows; row_idx > 0; --row_idx) {
        // Process the row for duplicates -- done in parallel
        row_processor(row_idx - 1);
    }
}

template <typename BaseBlock, size_t ThreadsX, size_t ThreadsY> 
void SubBlock<BaseBlock, ThreadsX, ThreadsY, devices::cpu>::determine_haplo_links()
{
    // Create a column processor to operate on the columns of the sub-block,
    // finding duplicate columns and determining the haplotype links
    Processor<sub_block_type, proc::col_dups_links, devices::cpu> col_processor(*this, _tree);
   
    // Start from the last column and go backwards 
    for (size_t col_idx = _cols; col_idx > 0; --col_idx) {
        // Set the number of elements for the nodes and their positions
        _tree.node(col_idx - 1).position() = col_idx - 1;
        _tree.node(col_idx - 1).elements() = _snp_info[col_idx - 1].length();
        
        // Check if the column is NIH, and set it if necessary
        if (_snp_info[col_idx - 1].type() == NIH) {
            _tree.node(col_idx - 1).type() = NIH;
            ++_num_nih;
        }
        
        // Process the column with the column processor to determine
        // duplicate columns and initialize the tree links and weights
        col_processor(col_idx - 1);
    }
}

}               // End namespace haplo
#endif          // PARAHAPLO_SUB_BLOCK_CPU_HPP


