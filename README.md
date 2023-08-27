# Simulation Study for Turning Circle: Ship Maneuver in Deep Water

This repository contains MATLAB code for simulating ship maneuvering during a turning circle in deep water. The provided MATLAB code allows you to study and visualize the turning behavior of a ship in a controlled simulation environment.

## Overview
Ship maneuverability is a critical aspect of maritime operations, and understanding how a ship responds during a turning circle is essential for safe navigation. This simulation study uses MATLAB to model the dynamics of a ship's turning circle in deep water, providing insights into factors such as turning radius, rate of turn, and stability during the maneuver.

### Key Features

1. **Hydrodynamic Coefficients & Forces:** This study incorporates hydrodynamic coefficients and forces, which play a significant role in the accuracy of the simulation. These coefficients are essential for modeling the complex interactions between the ship and the water.

2. **Equation of Motion and Time Integration:** The simulation is based on the equation of motion for ship maneuvering. The MATLAB code includes algorithms for time integration, which ensure that the simulation captures the ship's behavior accurately over time.

3. **New Rudder Calculation:** The simulation takes into account the ship's rudder dynamics, enabling the study of how rudder angles impact turning performance. This calculation enhances the realism of the simulation.

4. **IMO Criteria for Turning Circle Test:** To satisfy IMO criteria for this test, the advance (a) of the ship should not exceed 4.5 times its length (L), and the tactical diameter (d) should not exceed 5 times its length (L). These criteria are important for evaluating the ship's maneuverability within safe limits.

## Contents

- `turning_circle_simulation.m`: This MATLAB script implements the simulation of the ship's turning circle maneuver. It includes the necessary algorithms and equations to simulate the ship's dynamics and generate relevant data.

## How to Use

1. Open the `turning_circle_simulation.m` script and adjust the input parameters as needed. These parameters might include ship characteristics, initial conditions, and simulation duration.

2. Run the `turning_circle_simulation.m` script in MATLAB. This will perform the simulation and generate data. The visualize and analyze the simulation results. You can customize the analysis to focus on specific aspects of the turning circle maneuver.

## Contribution

Feel free to contribute to this repository by improving the simulation code, enhancing result analysis, or adding more visualization options. If you find any issues or have suggestions, please create an issue in the repository.

## License

This project is provided under the [MIT License](LICENSE), allowing you to use, modify, and distribute the code for both personal and commercial purposes.

## Acknowledgments

This simulation study is intended to enhance understanding of ship maneuverability and contribute to the field of maritime engineering. It serves as a starting point for exploring more complex ship dynamics and maneuvering scenarios.
