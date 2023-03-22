#include "Tissue.h"
#include "Cell.h"
#include "Initialization.h"
#include "SimulationBox.h"
// clang-format off
// The following includes must keep the order.
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <boost/numeric/odeint/integrate/integrate.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
// clang-format on
#include <cstddef>
#include <functional>
#include <oneapi/tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

Tissue::Tissue() { TissueInitializer::initializeTissue(this); }

Vector2F Tissue::computeInterCellForceOn1(const CellInd i_cell_1, const SubCellEnum i_subcell_1,
                                          const CellInd i_cell_2, const SubCellEnum i_subcell_2,
                                          const CellColonyMatrix &cols_of_cell_current_in) const {
    return CellBase::computeInterCellForceOn1(
        cols_of_cell_current_in.col(i_cell_1).segment<2>(CellBase::subcellBeginInd(i_subcell_1)),
        cols_of_cell_current_in.col(i_cell_2).segment<2>(CellBase::subcellBeginInd(i_subcell_2)));
}

Vector2F Tissue::computeInterCellForceOn1(const CellSubCellInd &i_cell_subcell_1,
                                          const CellSubCellInd &i_cell_subcell_2,
                                          const CellColonyMatrix &cols_of_cell_current_in) const {
    return computeInterCellForceOn1(i_cell_subcell_1.i_cell, i_cell_subcell_1.i_subcell,
                                    i_cell_subcell_2.i_cell, i_cell_subcell_2.i_subcell,
                                    cols_of_cell_current_in);
}

Vector2F
Tissue::computeInterCellForceOn1(const std::pair<CellSubCellInd, CellSubCellInd> &interacting_pair,
                                 const CellColonyMatrix &cols_of_cell_current_in) const {
    return computeInterCellForceOn1(interacting_pair.first, interacting_pair.second,
                                    cols_of_cell_current_in);
}

Matrix4XF
Tissue::computeInterCellForce_serial(const CellColonyMatrix &cols_of_cell_current_in) const {
    const CellInd tot_num_cells = states_of_cell_current.size();
    Matrix4XF inter_cell_force = Matrix4XF::Zero(Matrix4XF::RowsAtCompileTime, tot_num_cells);
    for (CellInd i_cell_1 = 0; i_cell_1 < tot_num_cells; ++i_cell_1) {
        for (CellInd i_cell_2 = i_cell_1 + 1; i_cell_2 < tot_num_cells; ++i_cell_2) {
            for (SubCellEnum i_subcell_1 : {SubCellEnum::front, SubCellEnum::back}) {
                for (SubCellEnum i_subcell_2 : {SubCellEnum::front, SubCellEnum::back}) {
                    const Vector2F inter_cell_force_on1 = computeInterCellForceOn1(
                        i_cell_1, i_subcell_1, i_cell_2, i_subcell_2, cols_of_cell_current_in);
                    inter_cell_force.col(i_cell_1).segment<2>(
                        CellBase::subcellBeginInd(i_subcell_1)) += inter_cell_force_on1;
                    inter_cell_force.col(i_cell_2).segment<2>(
                        CellBase::subcellBeginInd(i_subcell_2)) -= inter_cell_force_on1;
                }
            }
        }
    }
    return inter_cell_force;
}

Matrix4XF Tissue::computeInterCellForce(const CellColonyMatrix &cols_of_cell_current_in) const {
    const CellInd tot_num_cells = states_of_cell_current.size();
    std::vector<CellMutexes> all_cells_mutexes(tot_num_cells);
    Matrix4XF inter_cell_force = Matrix4XF::Zero(Matrix4XF::RowsAtCompileTime, tot_num_cells);

    auto compute_update_two_cells_inter_force = [&](const CellInd i_cell_1,
                                                    const CellInd i_cell_2) {
        for (SubCellEnum i_subcell_1 : {SubCellEnum::front, SubCellEnum::back}) {
            for (SubCellEnum i_subcell_2 : {SubCellEnum::front, SubCellEnum::back}) {
                const Vector2F inter_cell_force_on1 = computeInterCellForceOn1(
                    i_cell_1, i_subcell_1, i_cell_2, i_subcell_2, cols_of_cell_current_in);
                if (inter_cell_force_on1 != Vector2F::Zero()) {
                    {
                        SubCellMutex::scoped_lock a_sub_cell_scoped_lock(
                            all_cells_mutexes[i_cell_1].getASubCellMutex(i_subcell_1));
                        inter_cell_force.col(i_cell_1).segment<2>(
                            CellBase::subcellBeginInd(i_subcell_1)) += inter_cell_force_on1;
                    }
                    {
                        SubCellMutex::scoped_lock a_sub_cell_scoped_lock(
                            all_cells_mutexes[i_cell_2].getASubCellMutex(i_subcell_2));
                        inter_cell_force.col(i_cell_2).segment<2>(
                            CellBase::subcellBeginInd(i_subcell_2)) -= inter_cell_force_on1;
                    }
                }
            }
        }
    };

    tbb::parallel_for(tbb::blocked_range2d<CellInd>(0, tot_num_cells, 0, tot_num_cells),
                      [&](const tbb::blocked_range2d<CellInd> &cell_ind_range) {
                          for (CellInd i_cell_1 = cell_ind_range.rows().begin();
                               i_cell_1 < cell_ind_range.rows().end(); ++i_cell_1) {
                              for (CellInd i_cell_2 =
                                       std::max(i_cell_1 + 1, cell_ind_range.cols().begin());
                                   i_cell_2 < cell_ind_range.cols().end(); ++i_cell_2) {
                                  compute_update_two_cells_inter_force(i_cell_1, i_cell_2);
                              }
                          }
                      });
    return inter_cell_force;
}

// TODO hasn't been checked or used
Matrix4XF Tissue::computeExternalForce(const CellColonyMatrix &cols_of_cell_current_in,
                                       const FloatingType t) const {
    Eigen::Matrix<FloatingType, 1, Eigen::Dynamic> cells_center_y =
        (cols_of_cell_current_in.row(1) + cols_of_cell_current_in.row(3)) / 2.;
    const FloatingType y_cell_center_max = cells_center_y.maxCoeff();
    const FloatingType y_cell_center_min = cells_center_y.minCoeff();

    constexpr FloatingType artificial_traction_y_span = 50;
    constexpr FloatingType artificial_traction_magnitude =
        7; // referred to stall tension, which is 10

    const CellInd tot_num_cells = states_of_cell_current.size();
    Matrix4XF artificial_traction_force =
        Matrix4XF::Zero(Matrix4XF::RowsAtCompileTime, tot_num_cells);

    for (CellInd i_cell = 0; i_cell < tot_num_cells; ++i_cell) {
        if (cells_center_y(i_cell) > y_cell_center_max - artificial_traction_y_span) {
            // upper boundary pulling upward
            artificial_traction_force.col(i_cell).FRONT_POSITION_SEGMENT.y() =
                artificial_traction_force.col(i_cell).BACK_POSITION_SEGMENT.y() =
                    artificial_traction_magnitude;
        } else if (cells_center_y(i_cell) < y_cell_center_min + artificial_traction_y_span) {
            // lower boundary pulling downward
            artificial_traction_force.col(i_cell).FRONT_POSITION_SEGMENT.y() =
                artificial_traction_force.col(i_cell).BACK_POSITION_SEGMENT.y() =
                    -artificial_traction_magnitude;
        }
    }

    return artificial_traction_force;
}

void Tissue::computeODEdydt(const CellColonyMatrix &cols_of_cell_current_in,
                            CellColonyMatrix &cols_of_cell_dydt, const FloatingType /*t*/) const {
    const CellInd tot_num_cells = states_of_cell_current.size();
    Matrix4XF inter_cell_force = computeInterCellForce(cols_of_cell_current_in);

    //    code for running computeInterCellForce_serial and compare with
    //    parallel version
    /*
    std::cout << "inter_cell_force\n";
  //    std::cout << "parallel:\n" << inter_cell_force << '\n';
    Matrix4XF inter_cell_force_serial = computeInterCellForce_serial(cols_of_cell_current_in);
  //    std::cout << "serial:\n" << inter_cell_force_serial << '\n';
    FloatingType max_difference_force_parallel_serial =
            (inter_cell_force -
  inter_cell_force_serial).array().abs().maxCoeff(); if
  (max_difference_force_parallel_serial >= 1e-10) { std::cout <<
  "max_difference_force_parallel_serial = " <<
  max_difference_force_parallel_serial << '\n'; } else { std::cout << "OK\n";
    }
    assert(max_difference_force_parallel_serial < 1e-10);
  */

    tbb::parallel_for(static_cast<CellInd>(0), tot_num_cells, [&](const CellInd i_cell) {
        states_of_cell_current[i_cell]->computeODEdydt(cols_of_cell_current_in.col(i_cell),
                                                       inter_cell_force.col(i_cell),
                                                       cols_of_cell_dydt.col(i_cell));
    });
}

void Tissue::odeintAStep() {
    namespace odeint = boost::numeric::odeint;
    namespace plh = std::placeholders;
    auto start_range = cols_of_cell_current.data();
    const auto num_integration_eval = odeint::integrate<FloatingType>(
        std::bind(&Tissue::computeODEdydt, this, plh::_1, plh::_2, plh::_3), cols_of_cell_current,
        0., tissue_parameter.dt_full_step, tissue_parameter.dt_odeint);
    if (num_integration_eval >= 10) {
        std::cout << "odeint eval count = " << num_integration_eval << '\n';
    }
    //    may consider to use a fixed step size stepper to reduce integration
    //    function call overhead, and/or maybe tune dt
}

void Tissue::checkACell(const CellInd i_cell) const {
    const Eigen::Ref<const CellVector> cell_current_vector = cols_of_cell_current.col(i_cell);
    for (Eigen::Index i = 0; i < cell_current_vector.RowsAtCompileTime; ++i) {
        if (!std::isfinite(cell_current_vector(i))) {
            throw simulation_error("std::isfinite test failed for cell_current_vector of cell #" +
                                   std::to_string(i_cell));
        }
    }
}

void Tissue::updateACellState(
    const CellInd i_cell, std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
    CellColonyMatrix &cols_of_newborn_cell_current) {
    std::unique_ptr<CellBase> updated_cell_state = states_of_cell_current[i_cell]->updateStageState(
        cols_of_cell_current.col(i_cell), states_of_newborn_cell_current,
        cols_of_newborn_cell_current, tissue_parameter.dt_full_step);
    //    TODO support updateTypeState if needed
    if (updated_cell_state) {
        // The state of the current cell is changed and needs to be substituted back in.
        assert(typeid(*states_of_cell_current[i_cell]) != typeid(*updated_cell_state) &&
               "updated state is the same as before; wrong");
        states_of_cell_current[i_cell] = std::move(updated_cell_state);
    }
}

void Tissue::rewindWholeCellCoordinates() {
    tbb::parallel_for(static_cast<Eigen::Index>(0), cols_of_cell_current.cols(),
                      [&](const Eigen::Index i_cell) {
                          Eigen::Ref<Vector4F> a_cell_position =
                              cols_of_cell_current.col(i_cell).CELL_POSITION_SEGMENT;
                          CellBase::rewindWholeCellsPositionByBoundaryCondition(a_cell_position);
                      });
}

void Tissue::doAfterOdeint() {
    //    If a cell finishes dividing during this function call, that cell
    //    updates to be a child cell, and another child cell is appended to the
    //    sequences of_newborn_cell_current. The new cell has the same state as
    //    its sibling cell (the daughter cell occupying storage of the parent
    //    cell). i.e. it is already updated once appended. Cannot multi-thread
    //    because the containers are not thread-safe on cell division. i.e. It is
    //    not safe to have multiple threads adding cells simultaneously. Cannot
    //    multi-thread also because updateACellState may call global
    //    RandomHelper functions which are not thread-safe.

    const CellInd tot_num_cells_before_division = states_of_cell_current.size();
    assert(cols_of_cell_current.cols() == tot_num_cells_before_division);
    std::vector<std::unique_ptr<CellBase>> states_of_newborn_cell_current;
    CellColonyMatrix cols_of_newborn_cell_current;

    // update each cell state and give birth to new cells
    for (CellInd i_cell = 0; i_cell < tot_num_cells_before_division; ++i_cell) {
        states_of_cell_current[i_cell]->rewindSubCellsPositionByBoundaryCondition(
            cols_of_cell_current.col(i_cell).CELL_POSITION_SEGMENT);
        checkACell(i_cell);
        updateACellState(i_cell, states_of_newborn_cell_current, cols_of_newborn_cell_current);
    }

    const CellInd tot_num_newborn_cells = states_of_newborn_cell_current.size();
    assert(cols_of_newborn_cell_current.cols() == tot_num_newborn_cells);

    // merge newborn cells into current cell containers
    states_of_cell_current.insert(states_of_cell_current.end(),
                                  std::make_move_iterator(states_of_newborn_cell_current.begin()),
                                  std::make_move_iterator(states_of_newborn_cell_current.end()));
    cols_of_cell_current.conservativeResize(Eigen::NoChange,
                                            tot_num_cells_before_division + tot_num_newborn_cells);
    cols_of_cell_current.rightCols(tot_num_newborn_cells) = cols_of_newborn_cell_current;
    assert(states_of_cell_current.size() == cols_of_cell_current.cols());

    // rewind cells by boundary condition
    rewindWholeCellCoordinates();
}

void Tissue::doAStep() {
    assert(states_of_cell_current.size() == cols_of_cell_current.cols());
    odeintAStep();
    doAfterOdeint();
    assert(states_of_cell_current.size() == cols_of_cell_current.cols());
}

void Tissue::preRelax(std::size_t tot_num_relax_steps, TissueTrajectory &relaxation_trajectory,
                      CellCycleStageHistory &cells_stage_history) {
    assert(relaxation_trajectory.empty() && "relaxation_trajectory does not start out empty");

    for (CellInd i_cell = 0; i_cell < states_of_cell_current.size(); ++i_cell) {
        states_of_cell_current[i_cell] =
            std::make_unique<CellNonMotile>(cols_of_cell_current.col(i_cell));
    }

    // push_back initial config even if tot_num_relax_steps == 0
    relaxation_trajectory.push_back(cols_of_cell_current);
    cells_stage_history.push_back(collectCellsCycleStageRepr());

    for (std::size_t i_step = 1; i_step <= tot_num_relax_steps; ++i_step) {
        doAStep();
        relaxation_trajectory.push_back(cols_of_cell_current);
        cells_stage_history.push_back(collectCellsCycleStageRepr());
    }

    for (CellInd i_cell = 0; i_cell < states_of_cell_current.size(); ++i_cell) {
        states_of_cell_current[i_cell] =
            std::make_unique<CellContracting>(cols_of_cell_current.col(i_cell));
    }
}

std::string Tissue::collectCellsCycleStageRepr() const {
    std::string cells_stage_repr(states_of_cell_current.size(), ' ');
    for (CellInd i_cell = 0; i_cell < states_of_cell_current.size(); ++i_cell) {
        cells_stage_repr[i_cell] = states_of_cell_current[i_cell]->getCellCycleStageRepr();
    }
    assert(cells_stage_repr.find(' ') == std::string::npos);
    return cells_stage_repr;
}

void Tissue::doSimulation(std::size_t tot_num_simulation_steps,
                          TissueTrajectory &simulation_trajectory,
                          CellCycleStageHistory &cells_stage_history) {
    assert(simulation_trajectory.empty() && "simulation_trajectory does not start out empty");
    try {
        preRelax(ModelParameter::num_pre_relaxation_steps, simulation_trajectory,
                 cells_stage_history);
    } catch (const simulation_error &a_simulation_error) {
        std::cerr << a_simulation_error.what() << std::endl;
        std::cerr << "Simulation terminated early during pre-relaxing. "
                  << "Trajectory so far is returned." << std::endl;
    }

    for (std::size_t i_step = 1; i_step <= tot_num_simulation_steps; ++i_step) {
        try {
            doAStep();
            simulation_trajectory.push_back(cols_of_cell_current);
            cells_stage_history.push_back(collectCellsCycleStageRepr());
        } catch (const simulation_error &a_simulation_error) {
            std::cerr << a_simulation_error.what() << std::endl;
            std::cerr << "Simulation terminated early at step " << i_step
                      << ". Trajectory so far is returned." << std::endl;
            break;
        }
        /* catch (const tbb::captured_exception &a_tbb_captured_error) {
            std::cerr << a_tbb_captured_error.what() << std::endl;
            std::cerr << "Simulation terminated early at step " << i_step << ".
        Trajectory so far is returned."
                      << std::endl;
            break;
        }*/
    }
}
