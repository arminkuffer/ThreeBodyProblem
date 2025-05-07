
# Three-Body Simulation in C++

This project simulates the **three-body problem** in C++, where three spherical bodies interact under gravitational forces. The simulation includes two numerical integration methods: **Runge-Kutta** and **Leapfrog**. 

## Features
- Simulates the gravitational interactions between three spherical bodies in 2D.
- Includes two numerical methods: **Runge-Kutta** and **Leapfrog**.
- Visualizes the system using **SFML**.

## Libraries Used
- **SFML** for graphics rendering.
- **iostream**, **vector**, **cmath**, **deque** for basic functionality and math operations.

## Installation

### Prerequisites
- A C++ compiler (e.g., GCC, Clang)
- **SFML** library installed

### How To Install SFML
On **Linux** (e.g., Ubuntu):
```bash
sudo apt-get install libsfml-dev
```
**On Windows**: Follow the installation instructions on the SFML website.

### Compile the Code
To compile the project, run the following command:
```bash
g++ -o three-body-simulation main.cpp -lsfml-graphics -lsfml-window -lsfml-system
```
This will generate an executable file which you can then run.

### Switching Between Methods
By default, the simulation uses the Leapfrog method. If you'd like to experiment with the Runge-Kutta method instead, you can switch the methods by:

1. Commenting out the Leapfrog method code.

2. Uncommenting the Runge-Kutta method code.

Both methods can be found in the main() function of the code.

### Use Other Initial Conditions
By default, the simulation uses initial conditions that result in a figure-eight orbit. You can experiment with other initial positions 
in the code by changing the values of the initial positions (r1, r2, r3) and velocities (v1, v2, v3).


