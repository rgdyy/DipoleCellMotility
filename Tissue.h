#ifndef DIPOLECELL_TISSUE_H
#define DIPOLECELL_TISSUE_H

#include "CellBase.h"
#include "CommonTypes.h"
#include "ModelTissueParameter.h"
#include <Eigen/Core>
#include <memory>
#include <vector>

class Tissue {
    friend struct TissueInitializer;

  public:
    Tissue();

    static constexpr auto tissue_parameter = ModelParameter::tissue_parameter;

  private:
    std::vector<std::unique_ptr<CellBase>> states_of_cell_current;
    CellColonyMatrix cols_of_cell_current;

    Vector2F computeInterCellForceOn1(CellInd i_cell_1, SubCellEnum i_subcell_1, CellInd i_cell_2,
                                      SubCellEnum i_subcell_2,
                                      const CellColonyMatrix &cols_of_cell_current_in) const;

    Vector2F computeInterCellForceOn1(const CellSubCellInd &i_cell_subcell_1,
                                      const CellSubCellInd &i_cell_subcell_2,
                                      const CellColonyMatrix &cols_of_cell_current_in) const;

    Vector2F
    computeInterCellForceOn1(const std::pair<CellSubCellInd, CellSubCellInd> &interacting_pair,
                             const CellColonyMatrix &cols_of_cell_current_in) const;

    Matrix4XF computeInterCellForce_serial(const CellColonyMatrix &cols_of_cell_current_in) const;

    Matrix4XF computeInterCellForce(const CellColonyMatrix &cols_of_cell_current_in) const;

    Matrix4XF computeExternalForce(const CellColonyMatrix &cols_of_cell_current_in,
                                   const FloatingType t) const;

    // void computeODEdydt(const CellColonyMatrix &cols_of_cell_current_in,
    //                     CellColonyMatrix &cols_of_cell_dydt, FloatingType /*t*/) const;

    void computeODEdydt(const CellColonyMatrix &cols_of_cell_current_in,
                        CellColonyMatrix &cols_of_cell_dydt, FloatingType /*t*/) const;

    void odeintAStep();

    void checkACell(CellInd i_cell) const;

    void updateACellState(CellInd i_cell,
                          std::vector<std::unique_ptr<CellBase>> &states_of_newborn_cell_current,
                          CellColonyMatrix &cols_of_newborn_cell_current);

    void rewindWholeCellCoordinates();

    void doAfterOdeint();
    // rewinds according to boundary condition, checks each cell, updates cells' state, appends
    // newborn cells.

    void doAStep();

    void preRelax(std::size_t tot_num_relax_steps, TissueTrajectory &relaxation_trajectory,
                  CellCycleStageHistory &cells_stage_history);
    // begin with nonmotile state and relax, then set all to contracting state

    std::string collectCellsCycleStageRepr() const;

  public:
    void doSimulation(std::size_t tot_num_simulation_steps,
                                  TissueTrajectory &simulation_trajectory,
                                  CellCycleStageHistory &cells_stage_history);
};

#endif // DIPOLECELL_TISSUE_H
