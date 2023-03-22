#ifndef DIPOLECELL_INITIALCONDITION_H
#define DIPOLECELL_INITIALCONDITION_H

#include "CommonTypes.h"
#include "ModelCellParameter.h"

namespace InitialConditionParameter {
inline constexpr CellInd initial_num_rows = 30, initial_num_cols = 30;
inline constexpr FloatingType initial_cell_position_grid_size =
    ModelParameter::Shared::normal_cell_l0_full + 2;

inline constexpr bool enable_fourier_initial_variation = false;
inline constexpr FloatingType fourier_order_initial_variation = 4; // number of periods
inline constexpr FloatingType circular_frequency_initial_variation =
    2 * M_PI * fourier_order_initial_variation /
    initial_num_cols; // in the unit of (cell number)^-1
inline constexpr FloatingType amplitude_initial_variation =
    0.05; // in the unit of initial y-separation
} // namespace InitialConditionParameter

#endif // DIPOLECELL_INITIALCONDITION_H
