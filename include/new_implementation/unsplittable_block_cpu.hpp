/// @file   unsplittables_block_cpu.hpp
/// @brief  Header file for the unsplittable block class cpu implementation for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP
#define PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP

#include "unsplittable_block.hpp"

namespace haplo {
    
// Specialization for the CPU implementation of the unsplittable block
template <typename Block, size_t THI, size_t THJ>
class UnsplittableBlock<Block, THI, THJ, devices::cpu> : public Block {
public:
    // ------------------------------------------- ALIAS'S --------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
public:
    UnsplittableBlock(const Block& block);
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a pointer to the block which is a base class of this unsplittable block
    /// @return     A pointer the the block which is a base class of this unsplittable block
    // ------------------------------------------------------------------------------------------------------
    const Block* base_block() const { return static_cast<const Block*>(this); }
};

// -------------------------------------------- IMPLEMENTATIONS ---------------------------------------------

template <typename Block, size_t THI, size_t THJ>
UnsplittableBlock<Block, THI, THJ, devices::cpu>::UnsplittableBlock(const Block& block) 
: Block(block) 
{
}
    
}               // End namespace haplo
#endif          // PARAHAPLO_UNSPLITTABLE_BLOCK_CPU_HPP


