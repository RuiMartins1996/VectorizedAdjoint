#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iostream>

#include "../../lib/include/lib.hpp"

using namespace boost::numeric::odeint;
using namespace vectorizedadjoint;

// Define an observer to act as a logger to console
class Observer
{
  public:
    Observer() = default;

    template <typename State, typename Time>
    void operator()(const State &x, Time time)
    {
        std::cout << "t:" << time << "\t";
        std::cout << " x = [";
        for (auto &xi : x) {
            std::cout << xi << " ";
        }
        std::cout << " ]" << std::endl;
    }

  private:
    std::vector<double> times = {};
};

// System functor for the VanDerPol oscillator
class VanDerPol
{
  public:
    VanDerPol() = default;

    // Overload operator() to define the rhs of the system of ODEs
    template <typename T>
    void operator()(const std::vector<T> &x, std::vector<T> &dxdt, const std::vector<T> &mu, const T t) const
    {
        dxdt[0] = x[1];
        dxdt[1] = mu[0] * ((1.0 - x[0] * x[0]) * x[1] - x[0]);
    }
};

const double tol = 1e-5;

// Define the type of stepper to use to solve the system of ODEs
typedef odeint::runge_kutta_fehlberg78<std::vector<double>> stepper_type;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <tolerance>" << std::endl;
        return 1;
    }

    // Parse tolerance from command line arguments
    double tol = std::stod(argv[1]);

    std::vector<double> dxdt(2);
    std::vector<double> mu = {1e3};
    std::vector<double> x0 = {2.0, -2.0 / 3.0 + 10.0 / (81.0 * mu[0]) - 292.0 / (2187.0 * mu[0] * mu[0])};

    std::cout << "Initial conditions: x0 = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    // Time interval
    double ti = 0.0;
    double tf = 5e-1;
    // Time step (initial)
    double dt = 0.001;
    // Tolerances for the adaptive stepper
    double AbsTol = tol;
    double RelTol = tol;

    auto ctrlStepper = odeint::make_controlled<stepper_type>(AbsTol, RelTol);

    // System size
    int N = 2;
    // Number of parameters
    int Npar = 1;

    // Create a stepper object
    stepper_type stepper;
    // Create the system function instance
    auto vdp = VanDerPol();

    // Create a driver object
    Driver driver(N, N, Npar);

    // Forward pass
    size_t numSteps = runge_kutta(ctrlStepper, vdp, x0, mu, ti, tf, dt, driver);

    std::cout << "Number of steps: " << numSteps << std::endl;

    std::cout << "Solution: x = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    auto lambda = std::vector(N, std::vector<double>(N));

    lambda[0][0] = 1.0; // dx(t_f)/d x(t_f) = 1.0
    lambda[0][1] = 0.0; // dx(t_f)/d v(t_f) = 0.0
    lambda[1][0] = 0.0; // dv(t_f)/d x(t_f) = 0.0
    lambda[1][1] = 1.0; // dv(t_f)/d v(t_f) = 1.0

    auto muadj = std::vector(N, std::vector<double>(Npar));
    muadj[0][0] = 0.0; // dx(t_f)/d(mu) = 0
    muadj[1][0] = 0.0; // dv(t_f)/d(mu) = 0

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, muadj);

    // Inform driver of stepper butcher tableau, needed for reverse pass
    constructDriverButcherTableau(driver, stepper);
    // Record RHS function with automatic differentiation
    recordDriverRHSFunction(driver, vdp);

    adjointSolve(driver, mu);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Npar; j++) {
            std::cout << "mu[" << i << "][" << j << "] = " << muadj[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
