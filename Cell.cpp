#include "Cell.h"
#include "RandomHelper.h"
#include "SimulationBox.h"
#include <cassert>
#include <iostream>
#include <optional>

void CellContracting::initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) {
    resetElasticL0AndProbabilities(cell_current_vector,
                                   getCellParameter().cell_elastic_l0_contraction_full);
}

CellContracting::CellContracting(Eigen::Ref<CellVector> cell_current_vector) noexcept {
    initializeCellCurrentVector(cell_current_vector);
    //    OK here to call a virtual function since it is the final override.
}

CellCycleStageEnum CellContracting::getCellCycleStage() const {
    return CellCycleStageEnum::contracting;
}

char CellContracting::getCellCycleStageRepr() const {
    return static_cast<char>(CellCycleStageEnum::contracting);
}

void CellContracting::computeODEdydt_(
    const Eigen::Ref<const CellVector> &cell_current_vector,
    const Eigen::Ref<const Vector4F> &external_force_labxy, Eigen::Ref<CellVector> cell_dydt_vector,
    const FloatingType cell_elastic_coeff, const FloatingType free_contraction_speed,
    const FloatingType stall_tension, const FloatingType friction_coeff_front,
    const FloatingType friction_coeff_back, const FloatingType accum_rate_cil,
    const FloatingType confront_force_threshold_cil, const FloatingType decay_rate_cil,
    const FloatingType accum_rate_division, const FloatingType tension_threshold_division,
    const FloatingType decay_rate_division) {
    Matrix2F cols_of_cellxy_base;
    FloatingType cell_length;
    Vector4F external_force_cellxy;
    FloatingType cell_body_tension;
    Vector4F subcell_velocity_cellxy;
    computeSubCellVelocityLabxy(cell_current_vector, external_force_labxy, cols_of_cellxy_base,
                                cell_length, external_force_cellxy, cell_body_tension,
                                subcell_velocity_cellxy, cell_dydt_vector.CELL_POSITION_SEGMENT,
                                cell_elastic_coeff, friction_coeff_front, friction_coeff_back);

    refCellElasticL0(cell_dydt_vector) =
        -computeMotorContractionSpeed(cell_body_tension, stall_tension, free_contraction_speed);
    assert(getCellElasticL0(cell_dydt_vector) <= 0 &&
           "cell elastic l0 increasing during contraction; wrong");

    refUTurnProbability(cell_dydt_vector) =
        accum_rate_cil * std::max(0., -external_force_cellxy(0) - confront_force_threshold_cil) -
        decay_rate_cil;
    //    Note the sign: external_force_cellxy(0), external force on front
    //    subcell, is < 0 in collisions.

    refDivisionProbability(cell_dydt_vector) =
        accum_rate_division * (cell_body_tension - tension_threshold_division);
    // accum_rate_division * std::max(0., cell_body_tension - tension_threshold_division) -
    // decay_rate_division;
    // std::cout << cell_body_tension << '\n'; // TODO figure out why negative
}

void CellContracting::computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                                     const Eigen::Ref<const Vector4F> &external_force_labxy,
                                     Eigen::Ref<CellVector> cell_dydt_vector) const {
    const auto &param = getCellParameter();
    CellContracting::computeODEdydt_(
        cell_current_vector, external_force_labxy, cell_dydt_vector, param.cell_elastic_coeff,
        param.free_contraction_speed, param.stall_tension, param.friction_coeff_front,
        param.friction_coeff_back, param.accum_rate_cil, param.confront_force_threshold_cil,
        param.decay_rate_cil, param.accum_rate_division, param.tension_threshold_division,
        param.decay_rate_division);
}

std::optional<CellCycleStageEnum>
CellContracting::updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector,
                                   const FloatingType cell_elastic_l0_min,
                                   const bool cell_cil_enabled, const bool cell_division_enabled) {
    refDivisionProbability(cell_current_vector) =
        std::clamp(getDivisionProbability(cell_current_vector), 0., 1.);
    refUTurnProbability(cell_current_vector) =
        std::clamp(getUTurnProbability(cell_current_vector), 0., 1.);

    if (cell_division_enabled &&
        RandomHelper::bernoulli1D(getDivisionProbability(cell_current_vector))) {
        // If probabilistically decides to divide, then enter the dividing state.
        return CellCycleStageEnum::dividing;
    } else { // Not to divide
        if (getCellElasticL0(cell_current_vector) <= cell_elastic_l0_min) {
            // If it has contracted to the minimum allowed length, it should enter the protruding
            // state.
            return CellCycleStageEnum::protruding;
        } else { // Not to divide or protrude, continue to contract.
            if (cell_cil_enabled &&
                RandomHelper::bernoulli1D(getUTurnProbability(cell_current_vector))) {
                makeUTurn(cell_current_vector);
                std::cout << "A cell took a u-turn.\n";
            }
            // It is in contracting state right now and should remain contracting, so NO new state
            // to return.
            return std::nullopt;
        }
    }
}

std::optional<CellCycleStageEnum>
CellContracting::updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                  const FloatingType dt_full_step) {
    return CellContracting::updateCycleStage_(cell_current_vector,
                                              getCellParameter().cell_elastic_l0_contraction_min,
                                              ModelSwitch::cell_cil_uturn_enabled);
}

std::unique_ptr<CellBase> CellContracting::updateStageState(
    Eigen::Ref<CellVector> cell_current_vector,
    std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
    CellColonyMatrix &cols_of_newborn_cell_current, const FloatingType dt_full_step) {
    std::optional<CellCycleStageEnum> updated_stage =
        updateCycleStage(cell_current_vector, dt_full_step);
    if (updated_stage.has_value()) {
        switch (*updated_stage) {
        case CellCycleStageEnum::protruding:
            return std::make_unique<CellProtruding>(cell_current_vector);
        case CellCycleStageEnum::dividing:
            std::cout << "A cell wants to divide.\n";
            return std::make_unique<CellDividing>(cell_current_vector);
        default:
            assert(!"invalid updated cell cycle stage from CellContracting");
            throw std::logic_error("invalid updated cell cycle stage from CellContracting");
        }
    } else {
        return nullptr;
    }
}


void CellNonContracting::computeODEdydt_(const Eigen::Ref<const CellVector> &cell_current_vector,
                                         const Eigen::Ref<const Vector4F> &external_force_labxy,
                                         Eigen::Ref<CellVector> cell_dydt_vector,
                                         const FloatingType cell_elastic_coeff,
                                         const FloatingType cell_elastic_l0_in,
                                         const FloatingType friction_coeff_front,
                                         const FloatingType friction_coeff_back) {
    computeSubCellVelocityLabxy(cell_current_vector, external_force_labxy,
                                cell_dydt_vector.CELL_POSITION_SEGMENT, cell_elastic_coeff,
                                friction_coeff_front, friction_coeff_back);

    cell_dydt_vector.CELL_INTERNAL_SEGMENT.setZero();
    assert(getCellElasticL0(cell_current_vector) == cell_elastic_l0_in);
    assert(getUTurnProbability(cell_current_vector) == 0);
    assert(getDivisionProbability(cell_current_vector) == 0);
}

std::optional<CellCycleStageEnum>
CellExpanding::updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector,
                                 FloatingType cell_expand_target_length) {
    if (computeCellLengthSquared(cell_current_vector.CELL_POSITION_SEGMENT) >=
        cell_expand_target_length * cell_expand_target_length) {
        return CellCycleStageEnum::contracting;
    } else {
        return std::nullopt;
    }
}

void CellProtruding::initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) {
    resetElasticL0AndProbabilities(cell_current_vector,
                                   getCellParameter().cell_elastic_l0_protruding);
}

CellProtruding::CellProtruding(Eigen::Ref<CellVector> cell_current_vector) noexcept {
    initializeCellCurrentVector(cell_current_vector);
    //    OK here to call a virtual function since it is the final override.
}

CellCycleStageEnum CellProtruding::getCellCycleStage() const {
    return CellCycleStageEnum::protruding;
}

char CellProtruding::getCellCycleStageRepr() const {
    return static_cast<char>(CellCycleStageEnum::protruding);
}

void CellProtruding::computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                                    const Eigen::Ref<const Vector4F> &external_force_labxy,
                                    Eigen::Ref<CellVector> cell_dydt_vector) const {
    const auto &param = getCellParameter();
    CellNonContracting::computeODEdydt_(cell_current_vector, external_force_labxy, cell_dydt_vector,
                                        param.cell_elastic_coeff, param.cell_elastic_l0_protruding,
                                        param.friction_coeff_front,
                                        param.friction_coeff_back_protruding);
}

std::optional<CellCycleStageEnum>
CellProtruding::updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                 const FloatingType dt_full_step) {
    return CellExpanding::updateCycleStage_(cell_current_vector,
                                            getCellParameter().cell_elastic_l0_contraction_full);
}

std::unique_ptr<CellBase> CellProtruding::updateStageState(
    Eigen::Ref<CellVector> cell_current_vector,
    std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
    CellColonyMatrix &cols_of_newborn_cell_current, const FloatingType dt_full_step) {
    std::optional<CellCycleStageEnum> updated_stage =
        updateCycleStage(cell_current_vector, dt_full_step);
    if (updated_stage.has_value()) {
        switch (*updated_stage) {
        case CellCycleStageEnum::contracting:
            return std::make_unique<CellContracting>(cell_current_vector);
        default:
            assert(!"invalid updated cell cycle stage from CellProtruding");
            throw std::logic_error("invalid updated cell cycle stage from CellProtruding");
        }
    } else {
        return nullptr;
    }
}

void CellDividing::initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) {
    resetElasticL0AndProbabilities(cell_current_vector,
                                   getCellParameter().cell_elastic_l0_dividing);
}

CellDividing::CellDividing(Eigen::Ref<CellVector> cell_current_vector) noexcept {
    initializeCellCurrentVector(cell_current_vector);
    //    OK here to call a virtual function since it is the final override.
}


CellCycleStageEnum CellDividing::getCellCycleStage() const { return CellCycleStageEnum::dividing; }

char CellDividing::getCellCycleStageRepr() const {
    return static_cast<char>(CellCycleStageEnum::dividing);
}

void CellDividing::computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                                  const Eigen::Ref<const Vector4F> &external_force_labxy,
                                  Eigen::Ref<CellVector> cell_dydt_vector) const {
    const auto &param = getCellParameter();
    CellNonContracting::computeODEdydt_(cell_current_vector, external_force_labxy, cell_dydt_vector,
                                        param.cell_elastic_coeff, param.cell_elastic_l0_dividing,
                                        param.friction_coeff_back, param.friction_coeff_back);
    //    TODO Should friction coefficients be the same for both subcells in division?
}

std::optional<CellCycleStageEnum>
CellDividing::updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                               const FloatingType dt_full_step) {
    const auto &param = getCellParameter();
    return CellExpanding::updateCycleStage_(cell_current_vector,
                                            2 * param.cell_elastic_l0_contraction_full +
                                                param.divided_cells_separation);
}

std::unique_ptr<CellBase> CellDividing::updateStageState(
    Eigen::Ref<CellVector> cell_current_vector,
    std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
    CellColonyMatrix &cols_of_newborn_cell_current, const FloatingType dt_full_step) {
    std::optional<CellCycleStageEnum> updated_stage =
        updateCycleStage(cell_current_vector, dt_full_step);
    if (updated_stage.has_value()) {
        switch (*updated_stage) {
        case CellCycleStageEnum::contracting: {
            // final process to complete division, equivalent to
            // inserting two new subcells between the two
            // subcells of this cell

            assert(states_of_newborn_cell_current.size() == cols_of_newborn_cell_current.cols());

            // update child cells and subcells position
            CellVector newborn_cell_current_vector = cell_current_vector;
            slideSubCellAfterReplication(newborn_cell_current_vector.CELL_POSITION_SEGMENT,
                                         cell_current_vector.CELL_POSITION_SEGMENT);

            // append newborn cell vector and state
            cols_of_newborn_cell_current.conservativeResize(
                Eigen::NoChange, cols_of_newborn_cell_current.cols() + 1);
            cols_of_newborn_cell_current.rightCols<1>() = newborn_cell_current_vector;
            states_of_newborn_cell_current.push_back(
                std::make_unique<CellContracting>(newborn_cell_current_vector));

            assert(states_of_newborn_cell_current.size() == cols_of_newborn_cell_current.cols());
            std::cout << "A cell divided.\n";

            // return this cell's updated state
            return std::make_unique<CellContracting>(cell_current_vector);
        }
        default:
            assert(!"invalid updated cell cycle stage from CellDividing");
            throw std::logic_error("invalid updated cell cycle stage from CellDividing");
        }
    } else {
        return nullptr;
    }
}

void CellDividing::slideSubCellAfterReplication(
    Eigen::Ref<Vector4F> position_of_cell_to_slide_front,
    Eigen::Ref<Vector4F> position_of_cell_to_slide_back) const {
    assert(position_of_cell_to_slide_front == position_of_cell_to_slide_back &&
           "should be two copies of the position before sliding");

    Matrix2F cols_of_cellxy_base;
    const FloatingType cell_current_length =
        computeCellxyBaseGetCellLength(position_of_cell_to_slide_front, cols_of_cellxy_base);
    const Eigen::Ref<const Vector2F> cell_back_to_front_unit_vec = cols_of_cellxy_base.col(0);
    const FloatingType cell_target_length = getCellParameter().cell_elastic_l0_contraction_full;
    //     (cell_current_length - getCellParameter().divided_cells_separation) / 2;
    // assert(cell_target_length > 0 && "cell not long enough to cut into two cells; should have
    // been "
    //                                  "prevented from entering this step");

    position_of_cell_to_slide_front.FRONT_POSITION_SEGMENT =
        position_of_cell_to_slide_front.BACK_POSITION_SEGMENT +
        cell_target_length * cell_back_to_front_unit_vec;
    rewindSubCellsPositionByBoundaryCondition(position_of_cell_to_slide_front);

    position_of_cell_to_slide_back.BACK_POSITION_SEGMENT =
        position_of_cell_to_slide_back.FRONT_POSITION_SEGMENT -
        cell_target_length * cell_back_to_front_unit_vec;
    rewindSubCellsPositionByBoundaryCondition(position_of_cell_to_slide_back);
}

CellCycleStageEnum CellNonMotile::getCellCycleStage() const {
    return CellCycleStageEnum::nonmotile;
}

char CellNonMotile::getCellCycleStageRepr() const {
    return static_cast<char>(CellCycleStageEnum::nonmotile);
}

void CellNonMotile::computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                                   const Eigen::Ref<const Vector4F> &external_force_labxy,
                                   Eigen::Ref<CellVector> cell_dydt_vector) const {
    const auto &param = getCellParameter();
    CellNonContracting::computeODEdydt_(cell_current_vector, external_force_labxy, cell_dydt_vector,
                                        param.cell_elastic_coeff, param.cell_elastic_l0_nonmotile,
                                        param.friction_coeff_back, param.friction_coeff_back);
    //    Nonmotile cells have the same front/back friction coefficients.
}

std::optional<CellCycleStageEnum>
CellNonMotile::updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector) {
    return std::nullopt;
}

std::optional<CellCycleStageEnum>
CellNonMotile::updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                const FloatingType dt_full_step) {
    return std::nullopt;
}

void CellNonMotile::initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) {
    resetElasticL0AndProbabilities(cell_current_vector,
                                   getCellParameter().cell_elastic_l0_nonmotile);
}

CellNonMotile::CellNonMotile(Eigen::Ref<CellVector> cell_current_vector) noexcept {
    initializeCellCurrentVector(cell_current_vector);
    //    OK here to call a virtual function since it is the final override.
}

std::unique_ptr<CellBase> CellNonMotile::updateStageState(
    Eigen::Ref<CellVector> cell_current_vector,
    std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
    CellColonyMatrix &cols_of_newborn_cell_current, const FloatingType dt_full_step) {
    return nullptr;
}
