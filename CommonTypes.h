#ifndef DIPOLECELL_COMMONTYPES_H
#define DIPOLECELL_COMMONTYPES_H

#include <Eigen/Core>
#include <stdexcept>
#include <string>
#include <tbb/spin_mutex.h>
#include <vector>

constexpr int single_cell_vector_length = 7;

using CellInd = std::size_t;
using IntegerType = long;
// One should prefer class-specific index types for indexing, e.g. Eigen::Index
// for Eigen types.
using FloatingType = double;
using Vector2F = Eigen::Vector2d;
using Vector4F = Eigen::Vector4d;
using CellVector = Eigen::Matrix<FloatingType, single_cell_vector_length, 1>;
using VectorXF = Eigen::VectorXd;
using Array4F = Eigen::Array4d;
using ArrayXF = Eigen::ArrayXd;
using Matrix2F = Eigen::Matrix2d;
using Matrix2XF = Eigen::Matrix<FloatingType, 2, Eigen::Dynamic>;
using Matrix4XF = Eigen::Matrix<FloatingType, 4, Eigen::Dynamic>;
using CellColonyMatrix = Eigen::Matrix<FloatingType, single_cell_vector_length, Eigen::Dynamic>;
using MatrixXF = Eigen::MatrixXd;

using TissueTrajectory = std::vector<CellColonyMatrix>;
using CellCycleStageHistory = std::vector<std::string>;

// Definitions for meaningful segments of a single-cell vector
// see also CellBase::refSomething, CellBase::getSomething
#define FRONT_POSITION_SEGMENT head<2>()
#define BACK_POSITION_SEGMENT segment<2>(2)
#define CELL_POSITION_SEGMENT head<4>()
#define CELL_INTERNAL_SEGMENT tail<3>()
#define CELL_PROBABILITY_SEGMENT tail<2>()

enum class CellTypeEnum : char {
    normal = 'n',
};

enum class CellCycleStageEnum : char {
    nonmotile = 'n',
    contracting = 'c',
    protruding = 'p',
    dividing = 'd',
};

enum class SubCellEnum : char {
    front = '1',
    back = '2',
};

struct CellSubCellInd {
    CellInd i_cell;
    SubCellEnum i_subcell;
};

// Units: Length: microns; Time: min; Force: nN

struct CellParameter {
    FloatingType cell_elastic_coeff;
    FloatingType cell_elastic_l0_nonmotile;
    FloatingType cell_elastic_l0_contraction_full;
    FloatingType cell_elastic_l0_contraction_min;
    FloatingType cell_elastic_l0_protruding;
    FloatingType cell_elastic_l0_dividing;
    FloatingType free_contraction_speed;
    FloatingType stall_tension;
    FloatingType friction_coeff_front;
    FloatingType friction_coeff_back;
    FloatingType friction_coeff_back_protruding;
    FloatingType time_allowed_for_protrusion;
    FloatingType time_allowed_for_division;
    FloatingType divided_cells_separation;
    FloatingType accum_rate_cil;
    FloatingType confront_force_threshold_cil;
    FloatingType decay_rate_cil;
    FloatingType accum_rate_division;
    FloatingType tension_threshold_division;
    FloatingType decay_rate_division;
};

struct InterCellParameter {
    FloatingType lj_epsilon;
    FloatingType lj_sigma;
};

struct TissueParameter {
    FloatingType dt_full_step;
    FloatingType dt_odeint;
};

using SubCellMutex = tbb::spin_mutex;

struct CellMutexes {
    SubCellMutex front_mutex;
    SubCellMutex back_mutex;

    SubCellMutex &getASubCellMutex(SubCellEnum i_subcell) {
        switch (i_subcell) {
        case SubCellEnum::front:
            return front_mutex;
        case SubCellEnum::back:
            return back_mutex;
        default:
            assert(!"invalid subcell type");
            throw std::logic_error("invalid subcell type");
        }
    }
};

class simulation_error : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

#endif // DIPOLECELL_COMMONTYPES_H
