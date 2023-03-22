#include "CellBase.h"
#include <algorithm>
#include <cassert>
#include <iostream>

std::size_t CellBase::subcellBeginInd(const SubCellEnum i_subcell) {
    switch (i_subcell) {
    case SubCellEnum::front:
        return 0;
    case SubCellEnum::back:
        return 2;
    default:
        assert(!"invalid subcell type");
        throw std::logic_error("invalid subcell type");
    }
}

Vector2F CellBase::computeInterCellForceOn1(const Eigen::Ref<const Vector2F> &position_1,
                                            const Eigen::Ref<const Vector2F> &position_2) {
    return computeInterCellForceOn2(position_2, position_1);
}

Vector2F CellBase::computeInterCellForceOn2(const Eigen::Ref<const Vector2F> &position_1,
                                            const Eigen::Ref<const Vector2F> &position_2) {
    return computeInterCellForceOn2_<InterCellForce>(position_1, position_2);
}

void CellBase::rewindSubCellsPositionByBoundaryCondition(Eigen::Ref<Vector4F> cell_position) {
    SimulationBox::rewindPositionInPlace(cell_position.FRONT_POSITION_SEGMENT);
    SimulationBox::rewindPositionInPlace(cell_position.BACK_POSITION_SEGMENT);
}

void CellBase::rewindWholeCellsPositionByBoundaryCondition(Eigen::Ref<Vector4F> cell_position) {
    cell_position.FRONT_POSITION_SEGMENT =
        cell_position.BACK_POSITION_SEGMENT +
        SimulationBox::rewindDiffPosition(cell_position.FRONT_POSITION_SEGMENT -
                                          cell_position.BACK_POSITION_SEGMENT);
}

void CellBase::resetElasticL0AndProbabilities(Eigen::Ref<CellVector> cell_current_vector,
                                              const FloatingType cell_elastic_l0_in) {
    refCellElasticL0(cell_current_vector) = cell_elastic_l0_in;
    cell_current_vector.CELL_PROBABILITY_SEGMENT.setZero();
}

FloatingType CellBase::computeCellLengthSquared(const Eigen::Ref<const Vector4F> &cell_position) {
    return SimulationBox::rewindDiffPosition(cell_position.FRONT_POSITION_SEGMENT -
                                             cell_position.BACK_POSITION_SEGMENT)
        .squaredNorm();
}

FloatingType &CellBase::refCellElasticL0(Eigen::Ref<CellVector> cell_vector) {
    return cell_vector(4);
}

FloatingType CellBase::getCellElasticL0(const Eigen::Ref<const CellVector> &cell_vector) {
    return cell_vector(4);
}

FloatingType &CellBase::refUTurnProbability(Eigen::Ref<CellVector> cell_vector) {
    return cell_vector(5);
}

FloatingType CellBase::getUTurnProbability(const Eigen::Ref<const CellVector> &cell_vector) {
    return cell_vector(5);
}

FloatingType &CellBase::refDivisionProbability(Eigen::Ref<CellVector> cell_vector) {
    return cell_vector(6);
}

FloatingType CellBase::getDivisionProbability(const Eigen::Ref<const CellVector> &cell_vector) {
    return cell_vector(6);
}

FloatingType
CellBase::computeCellxyBaseGetCellLength(const Eigen::Ref<const Vector4F> &cell_position,
                                         Eigen::Ref<Matrix2F> cols_of_cellxy_base) {
    Vector2F cell_back_to_front = SimulationBox::rewindDiffPosition(
        cell_position.FRONT_POSITION_SEGMENT - cell_position.BACK_POSITION_SEGMENT);
    const FloatingType cell_length = cell_back_to_front.norm();
    Vector2F cell_base_vec_x = cell_back_to_front / cell_length;
    Vector2F cell_base_vec_y{-cell_base_vec_x.y(), cell_base_vec_x.x()};

    cols_of_cellxy_base << cell_base_vec_x, cell_base_vec_y;
    return cell_length;
}

void CellBase::castLabToCell(const Eigen::Ref<const Matrix2F> &cols_of_cellxy_base,
                             const Eigen::Ref<const Vector4F> &vec_quantity_labxy,
                             Eigen::Ref<Vector4F> vec_quantity_cellxy) {
    vec_quantity_cellxy << cols_of_cellxy_base.transpose() *
                               vec_quantity_labxy.FRONT_POSITION_SEGMENT,
        cols_of_cellxy_base.transpose() * vec_quantity_labxy.BACK_POSITION_SEGMENT;
}

void CellBase::castCellToLab(const Eigen::Ref<const Matrix2F> &cols_of_cellxy_base,
                             const Eigen::Ref<const Vector4F> &vec_quantity_cellxy,
                             Eigen::Ref<Vector4F> vec_quantity_labxy) {
    vec_quantity_labxy << cols_of_cellxy_base * vec_quantity_cellxy.FRONT_POSITION_SEGMENT,
        cols_of_cellxy_base * vec_quantity_cellxy.BACK_POSITION_SEGMENT;
}

FloatingType CellBase::computeCellBodyTension(FloatingType cell_length,
                                              const FloatingType cell_elastic_rest_len,
                                              const FloatingType cell_elastic_coeff) {
    // if (cell_length <= 0 || cell_length >= 2 * cell_elastic_rest_len) {
    //     // assert(!"cell body length out of handleable range");
    //     throw simulation_error("cell_length = " + std::to_string(cell_length) +
    //                            ", out of handleable range");
    // }

    cell_length =
        std::clamp(cell_length, 0.01 * cell_elastic_rest_len, 1.99 * cell_elastic_rest_len);
    FloatingType tension =
        -cell_elastic_coeff *
        (std::pow(cell_length, -1) + std::pow(cell_length - 2 * cell_elastic_rest_len, -1));
    if (cell_length < cell_elastic_rest_len) {
        assert(tension < 0);
    } else if (cell_length > cell_elastic_rest_len) {
        assert(tension > 0);
    }
    return tension;
}

FloatingType CellBase::computeMotorContractionSpeed(const FloatingType cell_body_tension,
                                                    const FloatingType stall_tension,
                                                    const FloatingType free_contraction_speed) {
    // TODO see the effect of molecular clutch
    // return free_contraction_speed;

    return cell_body_tension < stall_tension
               ? free_contraction_speed * (1 - cell_body_tension / stall_tension)
               : 0;
}

void CellBase::computeSubCellVelocityCellxy(const Eigen::Ref<const Vector4F> &external_force_cellxy,
                                            const FloatingType cell_body_tension,
                                            Eigen::Ref<Vector4F> subcell_velocity_cellxy,
                                            const FloatingType friction_coeff_front,
                                            const FloatingType friction_coeff_back) {
    //    Cell-xy has origin at back subcell and x-axis pointing from back to front.
    Vector4F intra_cell_force_cellxy{-cell_body_tension, 0, cell_body_tension, 0};
    Array4F sub_cell_friction_coeffs{friction_coeff_front, friction_coeff_front,
                                     friction_coeff_back, friction_coeff_back};
    subcell_velocity_cellxy =
        (external_force_cellxy + intra_cell_force_cellxy).array() / sub_cell_friction_coeffs;
}

void CellBase::computeSubCellVelocityLabxy(const Eigen::Ref<const CellVector> &cell_current_vector,
                                           const Eigen::Ref<const Vector4F> &external_force_labxy,
                                           Eigen::Ref<Vector4F> subcell_velocity_labxy,
                                           const FloatingType cell_elastic_coeff,
                                           const FloatingType friction_coeff_front,
                                           const FloatingType friction_coeff_back) {
    Matrix2F cols_of_cellxy_base;
    FloatingType cell_length;
    Vector4F external_force_cellxy;
    FloatingType cell_body_tension;
    Vector4F subcell_velocity_cellxy;
    computeSubCellVelocityLabxy(cell_current_vector, external_force_labxy, cols_of_cellxy_base,
                                cell_length, external_force_cellxy, cell_body_tension,
                                subcell_velocity_cellxy, subcell_velocity_labxy, cell_elastic_coeff,
                                friction_coeff_front, friction_coeff_back);
}

void CellBase::computeSubCellVelocityLabxy(
    const Eigen::Ref<const CellVector> &cell_current_vector,
    const Eigen::Ref<const Vector4F> &external_force_labxy,
    Eigen::Ref<Matrix2F> cols_of_cellxy_base, FloatingType &cell_length,
    Eigen::Ref<Vector4F> external_force_cellxy, FloatingType &cell_body_tension,
    Eigen::Ref<Vector4F> subcell_velocity_cellxy, Eigen::Ref<Vector4F> subcell_velocity_labxy,
    const FloatingType cell_elastic_coeff, const FloatingType friction_coeff_front,
    const FloatingType friction_coeff_back) {
    cell_length = computeCellxyBaseGetCellLength(
        cell_current_vector.CELL_POSITION_SEGMENT,
        cols_of_cellxy_base); // already considers periodic boundary condition
    castLabToCell(cols_of_cellxy_base, external_force_labxy, external_force_cellxy);
    cell_body_tension = computeCellBodyTension(
        cell_length, getCellElasticL0(cell_current_vector), cell_elastic_coeff);
    //    std::cout << "cell_body_tension = " << cell_body_tension << '\n';
    computeSubCellVelocityCellxy(external_force_cellxy, cell_body_tension,
                                           subcell_velocity_cellxy, friction_coeff_front,
                                           friction_coeff_back);
    castCellToLab(cols_of_cellxy_base, subcell_velocity_cellxy, subcell_velocity_labxy);
}

void CellBase::makeUTurn(Eigen::Ref<CellVector> cell_current_vector) {
    cell_current_vector.FRONT_POSITION_SEGMENT.swap(cell_current_vector.BACK_POSITION_SEGMENT);
    refUTurnProbability(cell_current_vector) = 0;
}
