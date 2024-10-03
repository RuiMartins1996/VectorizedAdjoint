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
    VectorizedAdjoint/aadc          The AADC library for automatic differentiation
    VectorizedAdjoint/lib           library implementation
    VectorizedAdjoint/examples      example programs
    VectorizedAdjoint/doc           user's guide

Dependencies:
============================================================
VectorizedAdjoint implementation requires boot::numeric::odeint, a library for ordinary differential equation (ODE) solvers within the Boost C++ Libraries. To install boost on Ubuntu, do:

    sudo apt update
    sudo apt install libboost-all-dev

Usage:
============================================================
For usage of VectorizedAdjoint, please refer to the user's guide, as well as the example programs containing drivers and makefiles which can be used as templates. If you are on Ubuntu, to run the user's guide in your default browser you can do:  
    cd doc/build/html
    xdg-open index.html
If this does not work you can try navigating to the indicated folder and drag and drop the index.html file on your preferred browser.


How to build the examples:
============================================================
    cd examples/ExampleFolder
    mkdir build
    cmake -S . -B ./build
    cd build
    make
    ./executable_name <command_line_arguments>

Link to preprint:
============================================================
Here is a [link](http://www.example.com) to the Arxiv preprint of the article regarding the development of this library and the introduction of SIMD vectorization in the discrete adjoint sensitivity analysis algorithm. 

Personal note:
============================================================
This is my first coding project. It is in a very early stage and suggestions and constructive criticism are very welcome. You can submit an issue on GitHub or email me at [rui.carlos.andrade.martins@gmail.com](mailto:rui.carlos.andrade.martins@gmail.com).