VectorizedAdjoint version 1.0.0
by Rui Martins
Released July 15, 2024

Description:
============================================================
VectorizedAdjoint is a C++ library that implements discrete adjoint sensitivity analysis for systems of parameter dependent ordinary differential equations solved using explicit adaptive Runge-Kutta methods. 
VectorizedAdjoint leverages vectorization and automatic differentiation to accelerate the computations.
VectorizedAdjoint has been tested with the gcc compiler.

Structure:
============================================================
VectorizedAdjoint contains the following directory structure:

    VectorizedAdjoint/README        general instructions
    VectorizedAdjoint/lib           library implementation
    VectorizedAdjoint/examples      example programs
    VectorizedAdjoint/doc           user's guide (TODO!)

Dependencies:
============================================================
VectorizedAdjoint implementation requires boot::numeric::odeint, a library for ordinary differential equation (ODE) solvers within the Boost C++ Libraries. To install boost on Ubuntu, do:

    sudo apt update
    sudo apt install libboost-all-dev

Usage:
============================================================
For usage of VectorizedAdjoint, please refer to the user's guide, as well as the example programs containing drivers and makefiles which can be used as templates. 


How to run the examples:
============================================================
    cd examples/VanDerPol
    mkdir build
    cmake -S . -B ./build
    cd build
    make
    ./vanderpol <tolerance>

    cd examples/GeneralizedLotkaVolterra:
    mkdir build
    cmake -S . -B ./build
    cd build
    make
    ./lotka <tolerance> <N>