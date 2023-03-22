# Dipole cell motility model

This repo contains C++/Python/Pybind11 code for simulating collective cell motility where each cell is a "dipole", i.e. represented by two active points. I open-sourced these codes because I have decided to instead focus on other research in my PhD study. As there would not likely be a scientific publication recently, I chose to publish the technical work.

## Model
You can think of the cells as diatomic molecules except that each "atom" is active, exerting active traction forces on the substrate which is assumed to have infinite mass. In addition to the traction force, there are also intra- and inter-cellular forces between subcell "atoms". Each cell goes through motility *cycles*, including many small steps of contraction, protrusion, and possibly cell division. 

## Code
- A object-oriented **state pattern** is adopted to implement state transition along motility cycles.
- Threading Building Blocks (TBB) is used for multi-threaded force calculation. Instead of using complicated lock-free molecular dynamics multi-threading strategy, I used a quick spin lock. Based on my experiments, the computation is fairly scalable up to at least 32 cores.
- Eigen and Boost-ODE are used to handle the math.
- Pybind11 is used to transfer C++ simulated trajectory to Python for postprocessing.

## How to build
Either install the conda env with `env.yaml`, or build the Docker image with `Dockerfile` (which used micromamba in Docker). With the conda env activated or inside the docker container, run `cmake -S . -B build -G Ninja && cmake --build build`; or use VSCode dev container extension.
