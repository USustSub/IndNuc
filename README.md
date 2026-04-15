# STM Model: Earthquake Nucleation in Depleting Gas Reservoirs

This repository contains a MATLAB-based seismo-thermomechanical (STM) model. It is designed to simulate the development of earthquakes on a pre-existing normal fault crosscutting a depleting gas reservoir. The reference model specifically investigates the Zeerijp region of the Groningen gas field.

* The main branch deals with a homogeneous underground, while the layeredModel branch deals with lithological layers.
* For an extensive code description, please go to `readme.pdf`.
* To cite the project, please refer to "Frictional healing and induced earthquakes on conventionally stable faults", Nature Communications (2023) by Meng Li, André R. Niemeijer, and Ylona van Dinther.

## Getting Started

To run the simulation:

1. Open MATLAB and navigate to the repository folder.
2. Run `Main-code.m`.
3. The script will automatically load parameters from `input.mat`, generate the staggered grid, and begin solving the partial differential equations (PDEs).

## Repository Structure

* `Main-code.m`: The core execution script. It handles grid generation, rate-and-state friction solving, adaptive time-stepping, and calculates stress/pressure changes.
* `generate-input.m` / `input.mat`: Configuration files containing all adjustable physical and numerical parameters.
* **Supporting Functions:** Various MATLAB files called by the main script to execute the simulation.

## Configuration & Inputs

All model parameters can be easily modified within the input files. Key adjustable parameters include:

* **Geometry:** Fault dip angle, model dimensions, and grid node counts.
* **Material Properties:** Densities (rock, fluid, gas), shear wave speed, Poisson ratio, and bulk modulus.
* **Friction (Rate-and-State):** Direct effect ($a$), evolution effect ($b$), characteristic slip distance ($L$), and dynamic slip weakening velocity.
* **Loading:** Far-field loading rate and pressure depleting rate.

## Outputs

The code generates multiple outputs to help visualize the fault's behavior over time:

* `output.txt`: A text file saved at every time step containing the elapsed model time, running time, min/max slip rates, and maximal displacement.
* **Data Checkpoints:** The script saves `.mat` files of the data at different time steps so users can reload data without rerunning the entire simulation.
* **Visuals:** The code automatically generates plots and a movie of the simulation, which can be saved at the end of the run.

## Model Physics

The simulation operates on a 2D plain strain model (depth 2000m to 4000m) using a fully staggered grid. The physics engine relies on:

* Rate-and-State Friction (Aging Law).
* Hooke's Law for calculating subsurface stress and strain.
* Lithostatic pressure gradients with hydrostatic fluid and gas pressures.
