#ifndef DIPOLECELL_SIMULATIONBOX_H
#define DIPOLECELL_SIMULATIONBOX_H

#include "CommonTypes.h"
#include <Eigen/Core>

class SimulationBox {
    friend class TissueInitializer;

    static constexpr bool periodic_boundary_x = false, periodic_boundary_y = false;
    inline static FloatingType x_size = 1, y_size = 1;

  public:
    static void rewindPositionInPlace(Eigen::Ref<Vector2F> position_coord);

    static Vector2F rewindDiffPosition(const Eigen::Ref<const Vector2F> &diff_position_coord);
};

#endif // DIPOLECELL_SIMULATIONBOX_H
