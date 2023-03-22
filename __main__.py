import time
import importlib
from PostProcessing import postProcess


def runSimulation(tot_num_steps, simulation_binary_dir="build"):
    print("Hello, Cell Motility Research World!")

    DipoleCell = importlib.import_module(f"{simulation_binary_dir}.DipoleCellPybind11")

    tissue = DipoleCell.Tissue()
    time_start = time.monotonic()
    simulation_trajectory = tissue.doSimulation(tot_num_steps)
    time_finish = time.monotonic()
    print(f"computation time elapsed: {time_finish - time_start} seconds")

    return simulation_trajectory


simulation_trajectory = runSimulation(20000)
postProcess(simulation_trajectory, save_trajectory_to_file=True)
