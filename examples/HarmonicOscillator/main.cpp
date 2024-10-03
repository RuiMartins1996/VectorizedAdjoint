#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iostream>

#include "../../lib/include/lib.hpp"

using namespace boost::numeric::odeint;

// System functor for the VanDerPol oscillator
class HarmonicOscillator
{
  private:
    double k = 1.0;

  public:
    HarmonicOscillator() = default;

    // Overload operator() to define the rhs of the system of ODEs
    template <typename T>
    void operator()(
        const std::vector<T> &r,
        std::vector<T> &drdt,
        const std::vector<T> &mu,
        const T t) const
    {
        drdt[0] = r[1];
        drdt[1] = -k * r[0] - mu[0] * r[1];
    }
};

// Define the type of stepper to use to solve the system of ODEs
typedef odeint::runge_kutta4<std::vector<double>> stepper_type;

int main(int argc, char *argv[])
{
    // System size
    int Nin = 2;
    // Number of cost functions
    int Nout = 3; //! DEBUG
    // Number of parameters
    int Npar = 1;

    // Create a driver object
    Driver driver(Nin, Nout, Npar);

    // Time interval and time step (initial)
    double ti = 0.0, tf = 10.0, dt = 0.01;

    std::vector<double> mu = {0.151};
    std::vector<double> r0 = {0.0, 1.0};

    // Create a stepper instance
    stepper_type stepper;
    // Create the system function instance
    auto hm = HarmonicOscillator();

    // Forward pass
    size_t numSteps = runge_kutta(stepper, hm, r0, mu, ti, tf, dt, driver);

    std::cout << "Number of steps: " << numSteps << std::endl;
    std::cout << "Solution: r = [" << r0[0] << ", " << r0[1] << "]" << std::endl;

    auto lambda = std::vector(Nout, std::vector<double>(Nin));

    // lambda[0][0] = 1.0 * r0[0]; // dE(t_f)dx(t_f) = k * x(t_f)
    // lambda[0][1] = 1.0 * r0[1]; // dE(t_f)dv(t_f) = m * v(t_f)

    //! DEBUG
    lambda[0][0] = 1.0;         // x(t_f)
    lambda[0][1] = 0.0;         //
    lambda[1][0] = 0.0;         //
    lambda[1][1] = 1.0;         // v(t_f)
    lambda[2][0] = 1.0 * r0[0]; // dE(t_f)dx(t_f) = k * x(t_f);
    lambda[2][1] = 1.0 * r0[1]; // dE(t_f)dv(t_f) = m * v(t_f);

    auto muadj = std::vector(Nout, std::vector<double>(Npar));
    muadj[0][0] = 0.0; // dE(t_f)d\mu = 0.0
    muadj[1][0] = 0.0;
    muadj[2][0] = 0.0;

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, muadj);

    // Inform driver of stepper butcher tableau, needed for reverse pass
    constructDriverButcherTableau(driver, stepper);
    // Record RHS function with automatic differentiation
    recordDriverRHSFunction(driver, hm);

    // Reverse pass to obtain the adjoints of the cost functions
    backpropagation::adjointSolve(driver, mu);

    std::cout << "dEdmr:" << lambda[2][0] << std::endl;
    std::cout << "dEdmv:" << lambda[2][1] << std::endl;
    std::cout << "dEdmu:" << muadj[0][0] << std::endl;

    for (int i = 0; i < Nout; i++) {
        for (int j = 0; j < Nin; j++) {
            std::cout << "lambda[" << i << "][" << j << "] = " << lambda[i][j] << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
