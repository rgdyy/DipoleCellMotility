import numpy as np
import matplotlib.pyplot as plt
import pathlib
from PostProcessing import load_trajectory_store_to_list, trajectory_list_to_3d_array

friction_coeff_front = (1.5,)
friction_coeff_back = (1.2,)
friction_coeff_back_protruding = (12,)

trajectory = trajectory_list_to_3d_array(
    load_trajectory_store_to_list(
        pathlib.Path.home() / "Downloads/simulation_trajectory_store.pickle"
    )
)

l0_trajectory = trajectory[:, 4, :]
# print(l0_trajectory.shape)


def bin_y_a_step(coord_a_step):
    cell_center_coordinates = (coord_a_step[0:2, :] + coord_a_step[2:4, :]) / 2.0
    cell_center_y = cell_center_coordinates[1, :].reshape(-1)
    num_bins = 20
    y_bins = np.histogram_bin_edges(cell_center_y, bins=num_bins)
    cells_y_bin_inds = np.digitize(cell_center_y, y_bins)
    return cells_y_bin_inds, num_bins


def average_l0_change_a_step(coord_a_step, l0_change):
    cells_y_bin_inds, num_bins = bin_y_a_step(coord_a_step)
    l0_change = np.abs(l0_change)
    average_l0_change = np.zeros(num_bins + 2)
    for i_bin in range(num_bins + 2):
        average_l0_change[i_bin] = np.mean(l0_change[cells_y_bin_inds == i_bin])

    return average_l0_change


i_step = 4000
plt.plot(
    average_l0_change_a_step(
        trajectory[i_step], l0_trajectory[i_step] - l0_trajectory[i_step - 1]
    )
)
plt.show()
