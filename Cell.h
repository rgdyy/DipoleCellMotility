#ifndef DIPOLECELL_CELL_H
#define DIPOLECELL_CELL_H

#include "CellBase.h"
#include "CommonTypes.h"
#include "ModelCellParameter.h"
#include <Eigen/Core>
#include <memory>
#include <optional>
#include <utility>

class CellContracting : public CellBase {
  private:
    void initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) final;

  public:
    explicit CellContracting(Eigen::Ref<CellVector> cell_current_vector) noexcept;

  public:
    CellCycleStageEnum getCellCycleStage() const final;

    char getCellCycleStageRepr() const final;

  private:
    static void
    computeODEdydt_(const Eigen::Ref<const CellVector> &cell_current_vector,
                    const Eigen::Ref<const Vector4F> &external_force_labxy,
                    Eigen::Ref<CellVector> cell_dydt_vector, FloatingType cell_elastic_coeff,
                    FloatingType free_contraction_speed, FloatingType stall_tension,
                    FloatingType friction_coeff_front, FloatingType friction_coeff_back,
                    FloatingType accum_rate_cil, FloatingType confront_force_threshold_cil,
                    FloatingType decay_rate_cil, FloatingType accum_rate_division,
                    FloatingType tension_threshold_division, FloatingType decay_rate_division);

  public:
    void computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                        const Eigen::Ref<const Vector4F> &external_force_labxy,
                        Eigen::Ref<CellVector> cell_dydt_vector) const final;

  private:
    static std::optional<CellCycleStageEnum>
    updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector, FloatingType cell_elastic_l0_min,
                      bool cell_cil_enabled = ModelSwitch::cell_cil_uturn_enabled,
                      bool cell_division_enabled = ModelSwitch::cell_division_enabled);
    // Only for CellContracting: this function also probabilistically does CIL.

  public:
    std::optional<CellCycleStageEnum> updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                                       FloatingType dt_full_step) final;

    std::unique_ptr<CellBase>
    updateStageState(Eigen::Ref<CellVector> cell_current_vector,
                     std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                     CellColonyMatrix &cols_of_newborn_cell_current,
                     FloatingType dt_full_step) final;
};

class CellNonContracting : public CellBase {
  protected:
    static void computeODEdydt_(const Eigen::Ref<const CellVector> &cell_current_vector,
                                const Eigen::Ref<const Vector4F> &external_force_labxy,
                                Eigen::Ref<CellVector> cell_dydt_vector,
                                FloatingType cell_elastic_coeff, FloatingType cell_elastic_l0_in,
                                FloatingType friction_coeff_front,
                                FloatingType friction_coeff_back);
};

class CellExpanding : public CellNonContracting {
  protected:
    static std::optional<CellCycleStageEnum>
    updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector,
                      FloatingType cell_expand_target_length);
    // returns contraction state if cell has expanded enough, i.e. has expanded to >=
    // cell_expand_target_length
};

class CellProtruding : public CellExpanding {
  private:
    void initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) final;

  public:
    explicit CellProtruding(Eigen::Ref<CellVector> cell_current_vector) noexcept;

  public:
    CellCycleStageEnum getCellCycleStage() const final;

    char getCellCycleStageRepr() const final;

  public:
    void computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                        const Eigen::Ref<const Vector4F> &external_force_labxy,
                        Eigen::Ref<CellVector> cell_dydt_vector) const final;

    std::optional<CellCycleStageEnum> updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                                       FloatingType dt_full_step) final;

    std::unique_ptr<CellBase>
    updateStageState(Eigen::Ref<CellVector> cell_current_vector,
                     std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                     CellColonyMatrix &cols_of_newborn_cell_current,
                     FloatingType dt_full_step) final;
};

class CellDividing : public CellExpanding {
  private:
    void initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) final;

  public:
    explicit CellDividing(Eigen::Ref<CellVector> cell_current_vector) noexcept;

  public:
    CellCycleStageEnum getCellCycleStage() const final;

    char getCellCycleStageRepr() const final;

  public:
    void computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                        const Eigen::Ref<const Vector4F> &external_force_labxy,
                        Eigen::Ref<CellVector> cell_dydt_vector) const final;

    std::optional<CellCycleStageEnum> updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                                       FloatingType dt_full_step) final;

    std::unique_ptr<CellBase>
    updateStageState(Eigen::Ref<CellVector> cell_current_vector,
                     std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                     CellColonyMatrix &cols_of_newborn_cell_current,
                     FloatingType dt_full_step) final;

  private:
    void slideSubCellAfterReplication(Eigen::Ref<Vector4F> position_of_cell_to_slide_front,
                                      Eigen::Ref<Vector4F> position_of_cell_to_slide_back) const;
    //    At the last step of cell division, two child cells initially have
    //    exactly the same state and position. This function adjusts cell length
    //    by sliding the subcells. It slides a subcell along the cell
    //    orientation while fixing the other subcell, so the cell gets to the
    //    target length.
};

class CellNonMotile : public CellNonContracting {
  private:
    void initializeCellCurrentVector(Eigen::Ref<CellVector> cell_current_vector) final;

  public:
    explicit CellNonMotile(Eigen::Ref<CellVector> cell_current_vector) noexcept;

  public:
    CellCycleStageEnum getCellCycleStage() const final;

    char getCellCycleStageRepr() const final;

  public:
    void computeODEdydt(const Eigen::Ref<const CellVector> &cell_current_vector,
                        const Eigen::Ref<const Vector4F> &external_force_labxy,
                        Eigen::Ref<CellVector> cell_dydt_vector) const final;

  private:
    static std::optional<CellCycleStageEnum>
    updateCycleStage_(Eigen::Ref<CellVector> cell_current_vector);

  public:
    std::optional<CellCycleStageEnum> updateCycleStage(Eigen::Ref<CellVector> cell_current_vector,
                                                       FloatingType dt_full_step) final;

    std::unique_ptr<CellBase>
    updateStageState(Eigen::Ref<CellVector> cell_current_vector,
                     std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                     CellColonyMatrix &cols_of_newborn_cell_current,
                     FloatingType dt_full_step) final;
};

#endif // DIPOLECELL_CELL_H
