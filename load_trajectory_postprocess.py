from PostProcessing import (
    postProcess,
    load_trajectory_store_to_list,
    load_cell_state_history_to_list,
)

postProcess(
    load_trajectory_store_to_list(),
    load_cell_state_history_to_list(),
    plot_initial_final=True,
    plot_interval=100,
    plot_animation=True,
)
