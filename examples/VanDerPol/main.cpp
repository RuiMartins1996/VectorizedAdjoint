#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../lib/include/backpropagation.hpp"
#include "../../lib/include/runge_kutta.hpp"

using namespace boost::numeric::odeint;

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

const double tol = 1e-5;

// Define the type of stepper to use to solve the system of ODEs
typedef std::vector<double> state_type;
typedef odeint::runge_kutta_fehlberg78<state_type> stepper_type;
// typedef odeint::runge_kutta_cash_karp54<state_type> stepper_type;
// typedef odeint::runge_kutta4_classic<state_type> stepper_type;
// typedef odeint::euler<state_type> stepper_type;

// System functor for the VanDerPol oscillator
template <typename T>
class VanDerPol
{
  public:
    VanDerPol() = default;

    // Overload operator() to define the rhs of the system of ODEs
    void operator()(const std::vector<T> &x, std::vector<T> &dxdt, std::vector<T> &mu, const T t) const
    {
        dxdt[0] = x[1];
        dxdt[1] = mu[0] * ((1.0 - x[0] * x[0]) * x[1] - x[0]);
    }
};

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <tolerance>" << std::endl;
        return 1;
    }

    // Parse tolerance from command line arguments
    double tol = std::stod(argv[1]);

    VanDerPol<double> vanderpol;

    std::vector<double> dxdt(2);
    std::vector<double> mu = {1e3};
    std::vector<double> x0 = {2.0, -2.0 / 3.0 + 10.0 / (81.0 * mu[0]) - 292.0 / (2187.0 * mu[0] * mu[0])};

    vanderpol(x0, dxdt, mu, 0.0);
    std::cout << "Initial conditions: x0 = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    // Time interval
    double ti = 0.0;
    double tf = 5e-1;
    // Time step (initial)
    double dt = 0.001;
    // Tolerances for the adaptive stepper
    double AbsTol = tol;
    double RelTol = tol;

    // System size
    int N = 2;
    // Number of parameters
    int Npar = 1;
    // Create a stepper object
    stepper_type stepper;

    void *driver_handle = new Driver(stepper, VanDerPol<idouble>(), N, N, Npar, x0);

    // Solve the system of ODEs
    // stepper_type rk4;
    size_t numSteps = runge_kutta(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        VanDerPol<double>(),
        x0, mu, ti, tf, dt, driver_handle,
        Observer());

    std::cout << "Number of steps: " << numSteps << std::endl;

    std::cout << "Solution: x = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    auto muMatrix = backpropagation::backpropagation(driver_handle, mu);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Npar; j++)
            std::cout << "mu[" << i << "][" << j << "] = " << muMatrix[i][j] << " ";
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}
