/// @file   unsplittables_block_cpu.hpp
/// @brief  Header file for the unsplittable block class cpu implementation for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP
#define PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP

#include "devices.hpp"
#include "equality_checker.hpp"
#include "processor_cpu.hpp"
#include "unsplittable_block.hpp"

#include <sstream>

namespace haplo {


// Specialization for the CPU implementation of the unsplittable block
template <typename BaseBlock, size_t THI, size_t THJ>
class UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu> : public BaseBlock {
public:
    // ------------------------------------------- ALIAS'S --------------------------------------------------`
    using ublock_type           = UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>;
    using binary_container      = BinaryVector<2>;              // Vector which uses 2 bits per element
    using atomic_vector         = std::vector<tbb::atomic<size_t>>;
    using node_container        = NodeContainer<devices::cpu>;
    using concurrent_umap       = typename BaseBlock::concurrent_umap;
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_t THREADS_I = THI;
    static constexpr size_t THREADS_J = THJ;
private:
    binary_container    _data;              //!< The data for the block
    size_t              _index;             //!< The index of the unsplittable block within the base block
    size_t              _start_idx;         //!< The index of the first column in the base block
    size_t              _end_idx;           //!< The index of the last column in the abse block
    size_t              _cols;              //!< The number of columns in the unsplittable block
    size_t              _rows;              //!< The number of rows in the unsplittable block
    
    template <typename FriendType, byte ProcessType, byte DeviceType>
    friend class Processor;
    
    concurrent_umap     _duplicate_rows;        
    concurrent_umap     _duplicate_cols;
    concurrent_umap     _monotone_cols;         //!< Monotone columns
    concurrent_umap     _nonih_cols;            //!< Non intrinsically heterozygous columns
    concurrent_umap     _row_multiplicities;    //!< The multiplicity of the rows
    concurrent_umap     _col_multiplicities;    //!< The multiplicity of the columns
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
    explicit UnsplittableBlock(const BaseBlock& block, const size_t index);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the element at position row_idx, col_idx
    /// @param[in]  row_idx     The index of the row of the element
    /// @param[in]  col_idx     The index of the column of the element
    // ------------------------------------------------------------------------------------------------------
    byte operator()(const size_t row_idx, const size_t col_idx) const 
    { 
        return _data.get(row_idx * _cols + col_idx);
    }
   
    // ---- DEBUGGING
    void print() const;
        
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a pointer to the block which is a base class of this unsplittable block
    /// @return     A pointer the the block which is a base class of this unsplittable block
    // ------------------------------------------------------------------------------------------------------
    const BaseBlock* base_block() const { return static_cast<const BaseBlock*>(this); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the data for the unsplittable block with the releavant data from the base block
    // ------------------------------------------------------------------------------------------------------
    void fill();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Find the duplicate columns, and initializes the nodes as they depend on the columns deing
    ///             duplicates (separate functions are not possible because of the performance implications
    // ------------------------------------------------------------------------------------------------------
    void find_duplicate_cols();    
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Find the duplicate rows
    // ------------------------------------------------------------------------------------------------------
    void find_duplicate_rows();

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines if a row is singul
    /// @param[in]  row_idx     The row to check for singularity
    /// @return     If the row is singular or not 
    // ------------------------------------------------------------------------------------------------------
    bool is_singular(const size_t row_idx) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets a row of data for the unsplittable block
    /// @param[in]  row_idx     The index of the row in the base block
    // ------------------------------------------------------------------------------------------------------
    void set_row_data(const size_t row_idx);
};

// ----------------------- DEBUGGING

template <typename BaseBlock, size_t THI, size_t THJ>
void UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::print() const 
{
      std::cout << "\n\n|--------------DATA------------|\n\n";
    
    for (size_t r = 0; r < _rows; ++r) {
        for (size_t c = 0; c < _cols; ++c) {
            if (_data.get(r * _cols +c) != 2 ) {
                std::cout << static_cast<unsigned>(_data.get(r * _cols + c)) << " ";
            } else {
                std::cout << "- ";
            }
        }
        std::cout << "\n";
    } 
    
    std::cout << "\n\n| ----------- ROW MPLICITIES -------------------|\n\n";
    
    for (size_t r = 0; r < _rows; ++r) {
        if (_row_multiplicities.find(r) != _row_multiplicities.end()) {
            std::cout << r << " : " << _row_multiplicities.at(r) << "\n";
        }
    }
    
    std::cout << "\n\n| ----------- ROW DUPS -------------------|\n\n";
    
    for (size_t r = 0; r < _rows; ++r) {
        if (_duplicate_rows.find(r) != _duplicate_rows.end()) {
            std::cout << r << " : " << _duplicate_rows.at(r) << "\n";
        }
    }
    std::cout << "\n\n| ----------- COL MPLICITIES -------------------|\n\n";
    
    for (size_t c = 0; c < _cols; ++c) {
        if (_col_multiplicities.find(c) != _col_multiplicities.end()) {
            std::cout << c << " : " << _col_multiplicities.at(c) << "\n";
        }
    }
    
    std::cout << "\n\n| ----------- COL DUPS -------------------|\n\n";
    
    for (size_t c = 0; c < _cols; ++c) {
        if (_duplicate_cols.find(c) != _duplicate_cols.end()) {
            std::cout << c << " : " << _duplicate_cols.at(c) << "\n";
        }
    }
    std::cout << "\n";            

}

// -------------------------------------------- IMPLEMENTATIONS ---------------------------------------------

template <typename BaseBlock, size_t THI, size_t THJ>
UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::UnsplittableBlock(const BaseBlock& block, 
                                                                        const size_t     index) 
: BaseBlock(block)                                                                      , 
  _data(0)                                                                              , 
  _index(index)                                                                         , 
  _start_idx(block.unsplittable_column(index))                                          ,
  _end_idx(block.unsplittable_column(index + 1))                                        ,  
  _cols(block.unsplittable_column(index + 1) - block.unsplittable_column(index) + 1)    ,
  _rows(0)
{
    std::ostringstream error_message;
    error_message   << "Index for unsplittable block past max index\n" 
                    << "\tindex : "     << index 
                    << "\tmax_index : " << block.num_unsplittable_blocks() - 1 << "\n"; 
    
   // Check that the index is valid 
   try{
       if (index > block.num_unsplittable_blocks() - 1) 
            throw std::out_of_range(error_message.str());
   } catch (const std::out_of_range& oor) { 
        std::cerr << "Out of Range error: " << oor.what() << '\n'; 
   }
   
   fill();
   find_duplicate_rows();
   find_duplicate_cols();
}

template <typename BaseBlock, size_t THI, size_t THJ> 
void UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::fill()
{
    // Go over each of the data rows and check for singularity
    for (size_t row_idx = 0; row_idx < BaseBlock::rows; ++row_idx) {
        if (!is_singular(row_idx)) {
            set_row_data(row_idx);
            ++_rows;
        }   
    }
   
}

template <typename BaseBlock, size_t THI, size_t THJ> 
void UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::find_duplicate_cols()
{
    node_container nodes(_cols);               // Nodes for the tree we are going to solve
    
    // Create a column processor to operate on the columns of this ublock
    Processor<ublock_type, proc::COL_DUP_LINKS, devices::cpu> col_processor(this);
    
    for (size_t col_idx = _cols; col_idx > 0; --col_idx) {
        // Set the node index
            
        // Process the column with the column processor to determine
        // the node connections and duplicate columns -- this is done in parallel
        col_processor(col_idx - 1, nodes);
        
        if (base_block()->is_monotone(col_idx + _start_idx - 1)) {                  // Monoton
            _monotone_cols[col_idx - 1] = 0;
        } else if (!base_block()->is_intrin_hetero(col_idx + _start_idx - 1)) {     // NIH
            _nonih_cols[col_idx - 1] = 0;
        }
    }
}

template <typename BaseBlock, size_t THI, size_t THJ>
void UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::find_duplicate_rows()
{
    // Create a processor for the rows to determine duplicates
    Processor<ublock_type, proc::ROW_DUP, devices::cpu> row_processor(this);
    
    // For each of the rows, from back to front
    for (size_t row_idx = _rows; row_idx > 0; --row_idx) {
        // Process the row for duplicates -- done in parallel
        row_processor(row_idx - 1);
    }
}

template <typename BaseBlock, size_t THI, size_t THJ> 
bool UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::is_singular(const size_t row_idx) const
{
    // Set the number of threads to use (we can use both dimensions threads)
    const size_t threads_x = (THI + THJ) < _cols ? (THI + THJ) : _cols;
    
    tbb::atomic<int> num_elements{0};
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, _cols, threads_x); 
                
                for (size_t it_x = 0; it_x < thread_iters_x  && num_elements <= 1; ++it_x) {
                    size_t col_idx  = _start_idx + it_x * threads_x + thread_idx;
                    
                    if (base_block()->operator()(row_idx, col_idx) <= 1) ++num_elements;
                }
            }
        }
    );
    return num_elements <= 1 ? true : false;
}

template <typename BaseBlock, size_t THI, size_t THJ> 
void UnsplittableBlock<BaseBlock, THI, THJ, devices::cpu>::set_row_data(const size_t row_idx)
{
    size_t current_data_size = _data.size();
    
    // Resize the data to add another row
    _data.resize(current_data_size + _cols);
    
    const size_t threads_x = (THI + THJ) < _cols ? (THI + THJ) : _cols;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, _cols, threads_x); 
                
                for (size_t it_x = 0; it_x < thread_iters_x ; ++it_x) {
                    size_t offset   = it_x * threads_x + thread_idx;
                    size_t col_idx  = _start_idx + offset;
                   
                    _data.set(current_data_size + offset, base_block()->operator()(row_idx, col_idx));
                }
            }
        }
    );
}

}               // End namespace haplo
#endif          // PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP


