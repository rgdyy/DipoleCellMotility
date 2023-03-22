#ifndef DIPOLECELL_CELLBASE_H
#define DIPOLECELL_CELLBASE_H

#include "CommonTypes.h"
#include "ModelCellParameter.h"
#include "SimulationBox.h"
#include <memory>
#include <optional>

class CellNormal;

class CellBase {
  public:
    virtual ~CellBase() = default;

    static std::size_t subcellBeginInd(SubCellEnum i_subcell);
    // the beginning index of subcell in a cell vector

    virtual CellCycleStageEnum getCellCycleStage() const = 0;

    virtual char getCellCycleStageRepr() const = 0;

  protected:
    static constexpr auto cell_parameter = ModelParameter::normal_cell_parameter;

  public:
    const CellParameter &getCellParameter() const { return cell_parameter; }

  protected:
    struct InterCellForce {
        static constexpr auto force_parameter = ModelParameter::inter_normal_normal_parameter;

        static FloatingType getCutoffDistance() {
            // could be used by computeInterCellForceOn2_, but only used by forceOn2 now
            return 2.5 * force_parameter.lj_sigma;
        }

        static Vector2F forceOn2(const Eigen::Ref<const Vector2F> &position_vector_1to2) {
            //  forceOn2 definition need not apply boundary rewind
            //  because it is done by computeInterCellForceOn2_

            if (position_vector_1to2.cwiseAbs().maxCoeff() >= getCutoffDistance()) {
                // If delta x or delta y >= cutoff_distance, the L2 distance is surely >=
                // cutoff_distance. For performance: less chance of computing r_squared.
                return Vector2F::Zero();
            } else {
                const FloatingType distance_1_2_squared = position_vector_1to2.squaredNorm();
                if (distance_1_2_squared >= std::pow(getCutoffDistance(), 2)) {
                    return Vector2F::Zero();
                } else {
                    FloatingType pow_6_frac_sigma_r = std::pow(
                        force_parameter.lj_sigma * force_parameter.lj_sigma / distance_1_2_squared,
                        3);
                    return 24 * force_parameter.lj_epsilon * position_vector_1to2 /
                           distance_1_2_squared * pow_6_frac_sigma_r * (2 * pow_6_frac_sigma_r - 1);
                }
            }
        }
    };

    template <class TInterCellForce>
    static Vector2F computeInterCellForceOn2_(const Eigen::Ref<const Vector2F> &position1,
                                              const Eigen::Ref<const Vector2F> &position2) {
        Vector2F position_vector_1to2 = SimulationBox::rewindDiffPosition(position2 - position1);
        // TInterCellForce::forceOn2 definition need not apply boundary rewind
        assert(position_vector_1to2 != Vector2F::Zero() &&
               "two subcells from different cells seem to overlap; maybe by mistake");
        return TInterCellForce::forceOn2(position_vector_1to2);
    }

  public:
    static Vector2F computeInterCellForceOn1(const Eigen::Ref<const Vector2F> &position_1,
                                             const Eigen::Ref<const Vector2F> &position_2);

    static Vector2F computeInterCellForceOn2(const Eigen::Ref<const Vector2F> &position_1,
                                             const Eigen::Ref<const Vector2F> &position_2);

  public:
    virtual void computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                                const Eigen::Ref<const Vector4F> &external_force_labxy,
                                Eigen::Ref<CellVector> cell_dydt_vector) const = 0;

    virtual std::optional<CellCycleStageEnum>
    updateCycleStage(Eigen::Ref<CellVector> cell_current_vector, FloatingType dt_full_step) = 0;

    virtual std::unique_ptr<CellBase>
    updateStageState(Eigen::Ref<CellVector> cell_current_vector,
                     std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                     CellColonyMatrix &cols_of_newborn_cell_current, FloatingType dt_full_step) = 0;
    /* returns the unique pointer to the updated state of this cell due to progressed cell motile
     * cycle stage. The function may modify the cell_current_vector including changing position and
     * other state variables. If the cell decides to remain in the same state, nullptr is returned.
     * If the cell decides to switch to another state, the unique pointer to the new cell state is
     * constructed and returned. If the cell has just finished dividing, i.e. switching to another
     * state from dividing, the newborn cell's vector and state are appended. The tissue will
     * substitute the state of the current cell and later concatenate the new cells to the cell
     * sequence.
     */

    static void rewindSubCellsPositionByBoundaryCondition(Eigen::Ref<Vector4F> cell_position);

    static void rewindWholeCellsPositionByBoundaryCondition(Eigen::Ref<Vector4F> cell_position);
    // such that each cell lie together, having length equal to the shortest distance between
    // subcells, even when crossing the periodic boundary.
    // should rewind immediately after odeint for a step

  protected:
    virtual void initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) = 0;

    static void resetElasticL0AndProbabilities(Eigen::Ref<CellVector> cell_current_vector,
                                               FloatingType cell_elastic_l0_in);

    static FloatingType computeCellLengthSquared(const Eigen::Ref<const Vector4F> &cell_position);
    //    There is no point in calculating this in cell xy-system. Thereby, the
    //    cell length is just x1 - x2.

    static FloatingType &refCellElasticL0(Eigen::Ref<CellVector> cell_vector);
    //    returns ref to the component of cell_vector for elastic rest length.
    //    It might be the length (LHS of the ODE) or its changing rate (RHS)

    static FloatingType getCellElasticL0(const Eigen::Ref<const CellVector> &cell_vector);
    //    Non-ref const getter version

    static FloatingType &refUTurnProbability(Eigen::Ref<CellVector> cell_vector);
    //    returns ref to the component of cell_vector for u-turn probability. It
    //    might be the probability (LHS of the ODE) or its changing rate (RHS)

    static FloatingType getUTurnProbability(const Eigen::Ref<const CellVector> &cell_vector);
    //    Non-ref const getter version

    static FloatingType &refDivisionProbability(Eigen::Ref<CellVector> cell_vector);
    //    returns ref to the component of cell_vector for cell division
    //    probability. It might be the probability (LHS of the ODE) or its
    //    changing rate (RHS)

    static FloatingType getDivisionProbability(const Eigen::Ref<const CellVector> &cell_vector);
    //    Non-ref const getter version

    static FloatingType
    computeCellxyBaseGetCellLength(const Eigen::Ref<const Vector4F> &cell_position,
                                   Eigen::Ref<Matrix2F> cols_of_cellxy_base);
    //    takes care of boundary condition; uses nearest images of two subcells
    //    in case of PBC. returns cell length, which is exactly the
    //    back-to-front subcell vector in cell-xy.

    static void castLabToCell(const Eigen::Ref<const Matrix2F> &cols_of_cellxy_base,
                              const Eigen::Ref<const Vector4F> &vec_quantity_labxy,
                              Eigen::Ref<Vector4F> vec_quantity_cellxy);

    static void castCellToLab(const Eigen::Ref<const Matrix2F> &cols_of_cellxy_base,
                              const Eigen::Ref<const Vector4F> &vec_quantity_cellxy,
                              Eigen::Ref<Vector4F> vec_quantity_labxy);

    static FloatingType computeCellBodyTension(FloatingType cell_length,
                                               FloatingType cell_elastic_rest_len,
                                               FloatingType cell_elastic_coeff);

    static FloatingType computeMotorContractionSpeed(FloatingType cell_body_tension,
                                                     FloatingType stall_tension,
                                                     FloatingType free_contraction_speed);
    //    always >=0; = 0 iff stalled

    static void computeSubCellVelocityCellxy(
        const Eigen::Ref<const Vector4F> &external_force_cellxy, FloatingType cell_body_tension,
        Eigen::Ref<Vector4F> subcell_velocity_cellxy, FloatingType friction_coeff_front,
        FloatingType friction_coeff_back);

    static void computeSubCellVelocityLabxy(const Eigen::Ref<const CellVector> &cell_current_vector,
                                            const Eigen::Ref<const Vector4F> &external_force_labxy,
                                            Eigen::Ref<Vector4F> subcell_velocity_labxy,
                                            FloatingType cell_elastic_coeff,
                                            FloatingType friction_coeff_front,
                                            FloatingType friction_coeff_back);

    static void computeSubCellVelocityLabxy(
        const Eigen::Ref<const CellVector> &cell_current_vector,
        const Eigen::Ref<const Vector4F> &external_force_labxy,
        Eigen::Ref<Matrix2F> cols_of_cellxy_base, FloatingType &cell_length,
        Eigen::Ref<Vector4F> external_force_cellxy, FloatingType &cell_body_tension,
        Eigen::Ref<Vector4F> subcell_velocity_cellxy, Eigen::Ref<Vector4F> subcell_velocity_labxy,
        FloatingType cell_elastic_coeff, FloatingType friction_coeff_front,
        FloatingType friction_coeff_back);
    //    A more verbose version with more in/out args.

    static void makeUTurn(Eigen::Ref<CellVector> cell_current_vector);
    //    swaps front/back position and reset U-Turn Probability to 0.
};

#endif // DIPOLECELL_CELLBASE_H
