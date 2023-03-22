import pathlib
import itertools
import pickle
import traceback
import multiprocessing
import numpy as np
import matplotlib as mpl
import matplotlib.figure as mplFigure
from mpl_toolkits.axes_grid1 import make_axes_locatable


def compute_a_step_cell_strain(simulation_trajectory, i_step):
    rows_of_cell_vector = simulation_trajectory[i_step].transpose()
    cell_length = np.linalg.norm(
        rows_of_cell_vector[:, [0, 1]] - rows_of_cell_vector[:, [2, 3]], axis=1
    )
    cell_l0 = rows_of_cell_vector[:, 4]
    cell_strain = (cell_length - cell_l0) / cell_l0
    return cell_strain


def mark_newborn_cells(simulation_trajectory, i_step, count_interval):
    """marks newborn cells since count_interval steps ago"""
    assert i_step >= 0
    is_newborn = np.zeros(simulation_trajectory[i_step].shape[1], int)
    if count_interval > 0 and i_step >= count_interval:
        i_prev_step = i_step - count_interval
        num_cells_prev = simulation_trajectory[i_prev_step].shape[1]
        num_cells_now = simulation_trajectory[i_step].shape[1]
        if num_cells_now > num_cells_prev:
            is_newborn[num_cells_prev:] = 1

    return is_newborn


class TrajectorySnapshot:
    def __init__(self, simulation_trajectory):
        self.simulation_trajectory = simulation_trajectory

    def plot_a_step(self, i_step, count_division_interval, savefig_filename=None):
        fig = mplFigure.Figure(figsize=(4.8 * 3.0, 6.4), constrained_layout=True)
        axs = fig.subplots(1, 3)

        self.plot_a_step_segments_uncolored(i_step, fig, axs[0])

        for i_ax in (1, 2):
            axs[i_ax].set_xlim(axs[0].get_xlim())
            axs[i_ax].set_ylim(axs[0].get_ylim())

        self.color_strain(i_step, fig, axs[1])

        self.color_division(i_step, count_division_interval, fig, axs[2])

        if savefig_filename is None:
            savefig_filename = f"plot_{i_step}"
        fig.savefig(f"{savefig_filename}.pdf")
        # plt.show()

    def plot_a_step_segments_uncolored(self, i_step, fig, ax):
        cell_position_1cell1row = self.simulation_trajectory[i_step].transpose()[:, :4]
        cell_position = cell_position_1cell1row.reshape((-1, 2, 2))

        cell_line_collection = ax.add_collection(
            mpl.collections.LineCollection(cell_position, color="tab:gray")
        )

        cell_back_point_collection = ax.scatter(
            cell_position[:, 1, 0], cell_position[:, 1, 1], s=9, c="tab:gray"
        )

        cell_front_point_collection = ax.scatter(
            cell_position[:, 0, 0], cell_position[:, 0, 1], s=9, c="tab:orange"
        )

        # ax.set_xlim((-25, 625))
        # ax.set_ylim((-250, 550))
        ax.set_aspect(1)

    def color_strain(self, i_step, fig, ax):
        cell_position_1cell1row = self.simulation_trajectory[i_step].transpose()[:, :4]
        cell_position = cell_position_1cell1row.reshape((-1, 2, 2))

        cell_line_collection = ax.add_collection(
            mpl.collections.LineCollection(cell_position)
        )

        # array is for color
        cell_line_collection.set_array(
            compute_a_step_cell_strain(self.simulation_trajectory, i_step)
        )
        cell_line_collection.set_cmap("Spectral")
        cell_line_collection.set_clim(vmin=-0.75, vmax=0.75)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="8%", pad=0.5)
        fig.colorbar(
            cell_line_collection, label="strain", cax=cax, orientation="horizontal"
        )

        ax.set_aspect(1)

    def color_division(self, i_step, count_interval, fig, ax):
        cell_position_1cell1row = self.simulation_trajectory[i_step].transpose()[:, :4]
        cell_position = cell_position_1cell1row.reshape((-1, 2, 2))

        cell_line_collection = ax.add_collection(
            mpl.collections.LineCollection(cell_position)
        )

        # array is for color
        cell_line_collection.set_array(
            mark_newborn_cells(self.simulation_trajectory, i_step, count_interval)
        )
        cell_line_collection.set_cmap("binary")
        cell_line_collection.set_clim(vmin=0, vmax=1)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="8%", pad=0.5)
        fig.colorbar(
            cell_line_collection, label="newborn", cax=cax, orientation="horizontal"
        )

        ax.set_aspect(1)


class TrajectoryAnimation:
    def __init__(self, simulation_trajectory, color_strain=False, color_division=False):
        self.simulation_trajectory = simulation_trajectory

        self.fig = mplFigure.Figure(figsize=(4.8, 6.4), constrained_layout=True)
        self.ax = self.fig.subplots()

        # if apply set_array, color is overriden, and colormap is used
        self.cell_line_collection = self.ax.add_collection(
            mpl.collections.LineCollection([], color="tab:gray")
        )
        self.cell_line_collection.set_animated(True)

        assert not (color_strain and color_division)
        self.color_strain = color_strain
        self.color_division = color_division
        if self.color_strain:
            self.cell_line_collection.set_cmap("Spectral")
            self.cell_line_collection.set_clim(vmin=-0.75, vmax=0.75)
            cblabel = "strain"
        if self.color_division:
            self.cell_line_collection.set_cmap("binary")
            self.cell_line_collection.set_clim(vmin=0, vmax=1)
            cblabel = "divided"

        if self.color_strain or self.color_division:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("top", size="8%", pad=0.5)
            self.fig.colorbar(
                self.cell_line_collection,
                label=cblabel,
                cax=cax,
                orientation="horizontal",
            )

        self.cell_front_point_collection = self.ax.scatter([], [], s=9, c="tab:orange")
        self.cell_front_point_collection.set_animated(True)

        self.cell_back_point_collection = self.ax.scatter([], [], s=9, c="tab:gray")
        self.cell_back_point_collection.set_animated(True)

        self.ax.set_xlim((-25, 625))
        self.ax.set_ylim((-250, 550))
        self.ax.set_aspect(1)

    def extract_a_step_cell_position(self, i_step):
        """extracts and reshapes i_step cell positions"""
        rows_of_cell_position = self.simulation_trajectory[i_step].transpose()[:, :4]
        return (
            rows_of_cell_position.reshape((-1, 2, 2)),
            rows_of_cell_position[:, :2],
            rows_of_cell_position[:, 2:],
        )

    def update_a_step_cells(self, i_step):
        """updates the plot for a frame of animation"""
        self.ax.set_title(f"time step {i_step}")
        (
            cell_position_grouped_by_cell,
            front_subcell_position,
            back_subcell_position,
        ) = self.extract_a_step_cell_position(i_step)
        self.cell_line_collection.set_segments(cell_position_grouped_by_cell)

        if self.color_strain:
            self.cell_line_collection.set_array(
                compute_a_step_cell_strain(self.simulation_trajectory, i_step)
            )
        if self.color_division:
            self.cell_line_collection.set_array(mark_newborn_cells(i_step))

        self.cell_front_point_collection.set_offsets(front_subcell_position)
        self.cell_back_point_collection.set_offsets(back_subcell_position)
        return (
            self.cell_line_collection,
            self.cell_front_point_collection,
            self.cell_back_point_collection,
        )

    def do_animation(self, movie_filename="trajectory"):
        ani = mpl.animation.FuncAnimation(
            self.fig,
            func=self.update_a_step_cells,
            frames=range(len(self.simulation_trajectory)),
            interval=50,
            blit=True,
        )
        ani.save(
            f"{movie_filename}"
            + ("_color_strain" if self.color_strain else "")
            + ("_color_division" if self.color_division else "")
            + ".mp4",
            dpi=200,
        )


def postProcess(
    simulation_trajectory,
    cells_stage_history,
    *,
    print_trajectory=False,
    save_trajectory_to_file=False,
    plot_initial_final=False,
    plot_interval=False,
    plot_animation=False,
):
    """plot_interval>=0, positive plot_interval implies plot_initial_final"""
    # Assert trajectory container's raw data is not copied
    # assert set(traj_a_time.flags['OWNDATA'] for traj_a_time in simulation_trajectory) == set((False,))

    exception_occurred = False

    if print_trajectory:
        try:
            print(simulation_trajectory, cells_stage_history)
        except Exception:
            exception_occurred = True
            traceback.print_exc()

    if save_trajectory_to_file:
        try:
            # np.savez_compressed(
            #     'simulation_trajectory_store', *simulation_trajectory)
            pickle.dump(
                simulation_trajectory, open("simulation_trajectory_store.pickle", "wb")
            )
            pickle.dump(
                cells_stage_history, open("cells_stage_history_store.pickle", "wb")
            )
        except Exception:
            exception_occurred = True
            traceback.print_exc()

    if plot_interval:
        plot_step_indices = range(0, len(simulation_trajectory), plot_interval)
        if plot_step_indices[-1] != len(simulation_trajectory) - 1:
            plot_step_indices = itertools.chain(
                plot_step_indices, (len(simulation_trajectory) - 1,)
            )
    elif plot_initial_final:
        plot_step_indices = (
            0,
            len(simulation_trajectory) - 1,
        )
    else:
        plot_step_indices = ()

    for i_step in plot_step_indices:
        try:
            trajectory_snapshot = TrajectorySnapshot(simulation_trajectory)
            trajectory_snapshot.plot_a_step(i_step, 1000)
            np.save(f"trajectory_{i_step}", simulation_trajectory[i_step])
        except Exception:
            exception_occurred = True
            traceback.print_exc()

    if plot_animation:
        try:
            trajectory_animation = TrajectoryAnimation(simulation_trajectory)
            # trajectory_animation_strain_colored = TrajectoryAnimation(
            #     simulation_trajectory, color_strain=True)
            # trajectory_animation_division_colored = TrajectoryAnimation(
            # simulation_trajectory, color_division=True)

            trajectory_animation.do_animation()

            # ### The following are for parallel processing of multiple movies, not used for now
            # with multiprocessing.Pool(2) as pool:
            #     ani0_proc = pool.apply_async(trajectory_animation.do_animation)
            #     # ani1_proc = pool.apply_async(
            #     #     trajectory_animation_strain_colored.do_animation)
            #     # ani2_proc = pool.apply_async(
            #     # trajectory_animation_division_colored.do_animation)

            #     ani0_proc.get()
            #     # ani1_proc.get()
            #     # ani2_proc.get()
            # ###
        except Exception:
            exception_occurred = True
            traceback.print_exc()

    assert not exception_occurred


def load_trajectory_store_to_list(store_file_path="simulation_trajectory_store.pickle"):
    store_file_path = pathlib.Path(store_file_path)
    if store_file_path.suffix == ".npz":
        trajectory_list = []
        with np.load(store_file_path) as trajectory_npz:
            for i_step in range(len(trajectory_npz.files)):
                trajectory_list.append(trajectory_npz[f"arr_{i_step}"])
        return trajectory_list
    elif store_file_path.suffix == ".pickle":
        with open(store_file_path, "rb") as trajectory_pickle:
            return pickle.load(trajectory_pickle)


def load_cell_state_history_to_list(
    history_store_path="cells_stage_history_store.pickle",
):
    history_store_path = pathlib.Path(history_store_path)
    with open(history_store_path, "rb") as history_pickle:
        return pickle.load(history_pickle)


def trajectory_list_to_3d_array(trajectory_list):
    """In the case without cell division, each step of trajectory should have the exact same shape.
    returns stacked trajectory indexed (time, coordinate, cell)"""
    return np.stack(trajectory_list, axis=0)
