#ifndef DIPOLECELL_MODELCELLPARAMETER_H
#define DIPOLECELL_MODELCELLPARAMETER_H

#include "CommonTypes.h"

namespace ModelSwitch {
inline constexpr bool cell_cil_uturn_enabled = false;
inline constexpr bool cell_division_enabled = false;
} // namespace ModelSwitch

namespace ModelParameter {
namespace Shared {
inline constexpr FloatingType normal_cell_l0_full = 5;
inline constexpr FloatingType normal_cell_stall_tension = 10;
} // namespace Shared

inline constexpr CellParameter normal_cell_parameter = {
    .cell_elastic_coeff = 10000,
    .cell_elastic_l0_nonmotile = Shared::normal_cell_l0_full,
    // .cell_elastic_l0_contraction_full = Shared::normal_cell_l0_full,
    .cell_elastic_l0_contraction_full = 0.65 * Shared::normal_cell_l0_full,
    // .cell_elastic_l0_contraction_min = 0.6 * Shared::normal_cell_l0_full,
    .cell_elastic_l0_contraction_min = 0.5 * Shared::normal_cell_l0_full,
    // .cell_elastic_l0_protruding = 1.5 * Shared::normal_cell_l0_full,
    .cell_elastic_l0_protruding = 1.25 * Shared::normal_cell_l0_full,
    .cell_elastic_l0_dividing = 2.5 * Shared::normal_cell_l0_full,
    .free_contraction_speed = 5,
    .stall_tension = Shared::normal_cell_stall_tension,
    .friction_coeff_front = 1.5,
    .friction_coeff_back = 1.2,
    .friction_coeff_back_protruding = 12,
    .time_allowed_for_protrusion = 0.1,
    .time_allowed_for_division = 0.5,
    .divided_cells_separation = 0.1 * Shared::normal_cell_l0_full,
    .accum_rate_cil = 0.01,
    .confront_force_threshold_cil = 0.8 * Shared::normal_cell_stall_tension,
    .decay_rate_cil = 0.01,
    .accum_rate_division = 1e-3,
    .tension_threshold_division = 0,
    // 0.8 * Shared::normal_cell_stall_tension,
    .decay_rate_division = 0.01,
};

inline constexpr InterCellParameter inter_normal_normal_parameter = {
    .lj_epsilon = 2,
    .lj_sigma = 15,
};
} // namespace ModelParameter

#endif // DIPOLECELL_MODELCELLPARAMETER_H
