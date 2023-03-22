#include "Initialization.h"
#include "Cell.h"
#include "CommonTypes.h"
#include "InitialCondition.h"
#include "ModelCellParameter.h"
#include "RandomHelper.h"
#include "SimulationBox.h"
#include "Tissue.h"
#include <cassert>
#include <cmath>
#include <iostream>

void TissueInitializer::initializeTissue(Tissue *tissue) {
    // initializeTissue1DX(tissue);
    initializeTissue2DTriangularLattice(tissue);
}

void TissueInitializer::initializeTissue2DTriangularLattice(Tissue *tissue) {
    using namespace InitialConditionParameter;

    // using initial_cell_class = CellNonMotile;

    constexpr CellInd initial_num_cells = initial_num_rows * initial_num_cols;
    tissue->cols_of_cell_current.resize(Eigen::NoChange, initial_num_cells);
    tissue->states_of_cell_current.resize(initial_num_cells);

    //    Triangular lattice
    Eigen::Matrix<FloatingType, 2, 3> cols_of_lattice_base_vector;
    cols_of_lattice_base_vector << Vector2F{1, 0}, Vector2F{0.5, 0.5 * std::sqrt(3)},
        Vector2F{-0.5, 0.5 * std::sqrt(3)};
    cols_of_lattice_base_vector *= initial_cell_position_grid_size;

    SimulationBox::x_size = initial_num_cols * cols_of_lattice_base_vector(0, 0);
    SimulationBox::y_size = initial_num_rows * cols_of_lattice_base_vector(1, 1);

    Vector2F leftmost_per_row{0, 0};
    // 0.5 * (cols_of_lattice_base_vector.col(0) + cols_of_lattice_base_vector.col(1));
    for (CellInd i_cell = 0, i_row = 0; i_row < initial_num_rows; ++i_row) {
        for (CellInd i_col = 0; i_col < initial_num_cols /* + i_row % 2*/; ++i_col) {
            // const Vector2F cell_midpoint_position =
            //     leftmost_per_row + i_col * cols_of_lattice_base_vector.col(0);
            // const FloatingType cell_initial_length =
            //     RandomHelper::realUniform1D(ModelParameter::Shared::_cell_l0_full - 0.5,
            //                                 ModelParameter::Shared::_cell_l0_full + 0.5);
            // const FloatingType cell_initial_orientation = RandomHelper::realUniform1D(0, 2 *
            // M_PI); FloatingType cell_initial_orientation; if (i_row == 0) {
            //     cell_initial_orientation = 1.5 * M_PI;
            // } else if (i_row == initial_num_rows - 1) {
            //     cell_initial_orientation = 0.5 * M_PI;
            // } else {
            //     cell_initial_orientation = RandomHelper::realUniform1D(0, 2 * M_PI);
            // }

            Vector2F cell_position_variation = Vector2F::Zero();
            // TODO figure out why /2. is ok but not /2
            if constexpr (enable_fourier_initial_variation) {
                cell_position_variation << 0,
                    (i_row - initial_num_rows / 2.) * amplitude_initial_variation *
                        cols_of_lattice_base_vector(1, 1) *
                        std::sin(circular_frequency_initial_variation * i_col);
            }

            const Vector2F cell_midpoint_position = leftmost_per_row +
                                                    i_col * cols_of_lattice_base_vector.col(0) +
                                                    cell_position_variation;

            const FloatingType cell_initial_length = ModelParameter::Shared::normal_cell_l0_full;

            const FloatingType cell_initial_orientation =
                i_row < initial_num_rows / 2
                    ? 1.5 * M_PI
                    : 0.5 * M_PI; // lower half pointing down; upper half up

            const Vector2F cell_midpoint_to_front =
                0.5 * cell_initial_length *
                Vector2F{std::cos(cell_initial_orientation), std::sin(cell_initial_orientation)};

            tissue->cols_of_cell_current.col(i_cell)
                << cell_midpoint_position + cell_midpoint_to_front,
                cell_midpoint_position - cell_midpoint_to_front, cell_initial_length, 0, 0;
            // tissue->states_of_cell_current[i_cell] =
            //     std::make_unique<initial_cell_class>(tissue->cols_of_cell_current.col(i_cell));

            ++i_cell;
        }
        leftmost_per_row += cols_of_lattice_base_vector.col(2 - i_row % 2);
        //        Odd i_row-th rows have one more cell than even i_row-th.
    }
}

void TissueInitializer::initializeTissue1DX(Tissue *tissue) {
    using initial_cell_class = CellContracting;

    constexpr CellInd initial_num_cells = 1;
    tissue->cols_of_cell_current.resize(Eigen::NoChange, initial_num_cells);
    tissue->states_of_cell_current.resize(initial_num_cells);

    constexpr FloatingType initial_cell_position_grid_size =
        0.5 * ModelParameter::inter_normal_normal_parameter.lj_sigma;
    const Vector2F unit_vector_along_1d_axis{1, 0};
    const Vector2F lattice_base_vector =
        initial_cell_position_grid_size * unit_vector_along_1d_axis;

    // These were used to test periodic boundary condition.
    // SimulationBox::x_size = initial_num_cells * initial_cell_position_grid_size;
    // SimulationBox::x_size = 2 * initial_cell_position_grid_size;

    Vector2F cell_midpoint_position{0, 0};
    for (CellInd i_cell = 0; i_cell < initial_num_cells; ++i_cell) {
        const FloatingType cell_initial_length = ModelParameter::Shared::normal_cell_l0_full;
        //        const FloatingType cell_initial_length =
        //        RandomHelper::realUniform1D(
        //                ModelParameter::Shared::_cell_l0_full - 0.5,
        //                ModelParameter::Shared::_cell_l0_full + 0.5);
        const FloatingType cell_initial_orientation = i_cell < initial_num_cells / 2 ? -1 : 1;
        const Vector2F cell_midpoint_to_front =
            0.5 * cell_initial_length * cell_initial_orientation * unit_vector_along_1d_axis;
        tissue->cols_of_cell_current.col(i_cell) << cell_midpoint_position + cell_midpoint_to_front,
            cell_midpoint_position - cell_midpoint_to_front, cell_initial_length, 0, 0;
        tissue->states_of_cell_current[i_cell] =
            std::make_unique<initial_cell_class>(tissue->cols_of_cell_current.col(i_cell));

        cell_midpoint_position += lattice_base_vector;
    }
}
