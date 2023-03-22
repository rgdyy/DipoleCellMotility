#include "SimulationBox.h"
#include "CommonTypes.h"
#include "InitialCondition.h"
#include <cmath>
#include <functional>

// TODO check periodic boundary condition

void SimulationBox::rewindPositionInPlace(Eigen::Ref<Vector2F> position_coord) {
    //    https://en.wikipedia.org/wiki/Periodic_boundary_conditions#(A)_Restrict_particle_coordinates_to_the_simulation_box
    if constexpr (periodic_boundary_x) {
        position_coord.x() -= std::floor(position_coord.x() / x_size) * x_size;
    }
    if constexpr (periodic_boundary_y) {
        position_coord.y() -= std::floor(position_coord.y() / y_size) * y_size;
    }
}

Vector2F SimulationBox::rewindDiffPosition(const Eigen::Ref<const Vector2F> &diff_position_coord) {
    //    https://en.wikipedia.org/wiki/Periodic_boundary_conditions#(A)_Restrict_particle_coordinates_to_the_simulation_box
    Vector2F shifted_diff_position = diff_position_coord;
    if constexpr (periodic_boundary_x) {
        shifted_diff_position.x() -= std::nearbyint(shifted_diff_position.x() / x_size) * x_size;
    }
    if constexpr (periodic_boundary_y) {
        shifted_diff_position.y() -= std::nearbyint(shifted_diff_position.y() / y_size) * y_size;
    }
    return shifted_diff_position;
}
