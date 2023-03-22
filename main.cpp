#include "PostProcessingPybind11.h"
#include "Tissue.h"
#include <chrono>
#include <iostream>
#include <vector>

void runSimulation(std::size_t tot_num_steps, TissueTrajectory &simulation_trajectory,
                   CellCycleStageHistory &cells_stage_history) {
    std::cout << "Hello, Cell Motility Research World!" << std::endl;

    // auto tissue = Tissue();
    Tissue tissue{};
    const auto time_start = std::chrono::steady_clock::now();
    tissue.doSimulation(tot_num_steps, simulation_trajectory, cells_stage_history);
    const auto time_finish = std::chrono::steady_clock::now();
    const std::chrono::duration<double> computation_time_elapsed = time_finish - time_start;
    std::cout << "computation_time_elapsed = " << computation_time_elapsed.count() << " seconds"
              << std::endl;

    //    constexpr auto sep = "\n----------------------------------------\n";
    //    for (std::size_t i_step = 0; auto &position_a_time :
    //    simulation_trajectory) {
    //        std::cout << "Step " << i_step++ << sep << position_a_time << sep;
    //    }
}

int main() {
    TissueTrajectory simulation_trajectory;
    CellCycleStageHistory cells_stage_history;
    runSimulation(5000, simulation_trajectory, cells_stage_history);
    invokePostProcess(std::move(simulation_trajectory), std::move(cells_stage_history));
    return 0;
}
