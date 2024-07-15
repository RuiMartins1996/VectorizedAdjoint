# VectorizedAdjoint
This is a C++ library that implements discrete adjoint sensitivity analysis for systems of ODEs solved using explicit adaptive Runge-Kutta methods. This library leverages vectorization and automatic differentiation to accelerate the computations.


# How to run the examples
VanDerPol:
mkdir build
cmake -S . -B ./build
cd build
make
./vanderpol <tolerance>

GeneralizedLotkaVolterra:
mkdir build
cmake -S . -B ./build
cd build
make
./lotka <tolerance> <N>