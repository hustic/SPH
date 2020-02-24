# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Compilation/Installation Guide

*To be written*

### User instructions

*To be written*

### Documentation

*To be written*

### Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
python run_tests.py
```
