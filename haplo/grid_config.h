// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo grid configuration for cuda grids
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRID_CONFIG_H
#define PARHAPLO_GRID_CONFIG_H

namespace haplo {
namespace grid  {
    
static constexpr size_t block_size_x = 1024;
static constexpr size_t block_size_y = 1024;
static constexpr size_t block_size_z = 1024;


}               // End namespace grid
}               // End namespace haplo
#endif          // PARAHAPLO_GRID_CONFIG_H
