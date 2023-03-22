#include "PostProcessingPybind11.h"
#include "CommonTypes.h"
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <utility>

void invokePostProcess(TissueTrajectory simulation_trajectory, CellCycleStageHistory cells_stage_history) {
    pybind11::scoped_interpreter guard{};
    using namespace pybind11::literals;

    pybind11::object postProcess = pybind11::module_::import("PostProcessing").attr("postProcess");
    postProcess(std::move(simulation_trajectory), std::move(cells_stage_history), "save_trajectory_to_file"_a = true);
}
