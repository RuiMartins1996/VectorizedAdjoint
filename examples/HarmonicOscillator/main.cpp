#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iostream>

#include "../../lib/include/lib.hpp"

using namespace boost::numeric::odeint;

// System functor for the VanDerPol oscillator
class HarmonicOscillator
{
  public:
    HarmonicOscillator() = default;

    // Overload operator() to define the rhs of the system of ODEs
    template <typename T>
    void operator()(
        const std::vector<T> &x,
        std::vector<T> &dxdt,
        const std::vector<T> &mu,
        const T t) const
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - mu[0] * x[1];
    }
};

// Define the type of stepper to use to solve the system of ODEs
typedef odeint::runge_kutta4<std::vector<double>> stepper_type;

int main(int argc, char *argv[])
{
    std::vector<double> mu = {0.15};
    std::vector<double> x0 = {0.0, 1.0};

    // Time interval
    double ti = 0.0;
    double tf = 10.0;
    // Time step (initial)
    double dt = 0.01;

    // Create a stepper object
    stepper_type stepper;
    // Create the system function
    auto hm = HarmonicOscillator();

    // System size
    int N = 2;
    // Number of parameters
    int Npar = 1;
    // Create a driver object
    Driver driver(N, N, Npar);

    // Forward pass
    size_t numSteps = runge_kutta(stepper, hm, x0, mu, ti, tf, dt, driver);

    std::cout << "Number of steps: " << numSteps << std::endl;
    std::cout << "Solution: x = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    // Inform driver of stepper butcher tableau, needed for reverse pass
    constructDriverButcherTableau(driver, stepper);
    // Record RHS function with automatic differentiation
    recordDriverRHSFunction(driver, hm);

    // Compute the ODE sensitivity matrix
    auto muMatrix = backpropagation::computeSensitivityMatrix(driver, mu);

    // Print sensitivity matrix
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Npar; j++)
            std::cout << "mu[" << i << "][" << j << "] = " << muMatrix[i][j] << " ";
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
