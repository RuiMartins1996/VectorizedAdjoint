#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <boost/numeric/odeint.hpp>

#include "../../../lib/include/lib.hpp"

using namespace boost::numeric::odeint;
using namespace vectorizedadjoint;

const double tol = 1e-5;

// Define the type of stepper to use to solve the system of ODEs
typedef std::vector<double> state_type;
typedef runge_kutta_cash_karp54<state_type> stepper_type;

// System function for the VanDerPol oscillator
class VanDerPol
{
public:
    VanDerPol() = default;

    // Overload operator() to define the rhs of the system of ODEs
    template <typename T>
    void operator()(const std::vector<T> &x, std::vector<T> &dxdt,
                    std::vector<T> &eps, const T t) const
    {
        dxdt[0] = x[1];
        dxdt[1] = ((1.0 - x[0] * x[0]) * x[1] - x[0])/eps[0];
    }
};

int main(int argc, char* argv[])
{
        if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <tolerance>" << std::endl;
        return 1;
    }

    // Parse tolerance from command line arguments
    double tol = std::stod(argv[1]);

    std::vector<double> dxdt(2);
    std::vector<double> eps = {0.001};
    std::vector<double> x0 = {2.0, -2.0 / 3.0 + 10.0 / (81.0 / eps[0]) - 292.0 / (2187.0 / (eps[0] * eps[0]))};

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
    size_t numSteps = runge_kutta(ctrlStepper, vdp, x0, eps, ti, tf, dt, driver);

    std::cout << "Number of steps: " << numSteps << std::endl;

    std::cout << "Solution: x = [" << x0[0] << ", " << x0[1] << "]" << std::endl;

    auto lambda = std::vector(N, std::vector<double>(N));

    lambda[0][0] = 1.0; // dx(t_f)/d x(t_f) = 1.0
    lambda[0][1] = 0.0; // dx(t_f)/d v(t_f) = 0.0
    lambda[1][0] = 0.0; // dv(t_f)/d x(t_f) = 0.0
    lambda[1][1] = 1.0; // dv(t_f)/d v(t_f) = 1.0

    auto epsadj = std::vector(N, std::vector<double>(Npar));
    epsadj[0][0] = 0.0; // dx(t_f)/d(eps) = 0
    epsadj[1][0] = 0.0; // dv(t_f)/d(eps) = 0

    // Set derivatives of cost functions w.r.t ODE solution and w.r.t. parameters
    setCostGradients(driver, lambda, epsadj);

    // Inform driver of stepper butcher tableau, needed for reverse pass
    constructDriverButcherTableau(driver, stepper);
    // Record RHS function with automatic differentiation
    recordDriverRHSFunction(driver, vdp);

    adjointSolve(driver, eps);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Npar; j++) {
            std::cout << "eps[" << i << "][" << j << "] = " << epsadj[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Print to CSV file
    std::ofstream outfile("output.csv");
    if (outfile.is_open()) {
        // Write header
        outfile << "x0, x1, epsAdjoints0, epsAdjoints1" << std::endl;
        
        outfile << std::scientific << std::setprecision(16);
        // Write values
        outfile << x0[0] << ", " << x0[1] << ", " << epsadj[0][0] << ", " << epsadj[1][0] << std::endl;

        outfile.close();
        std::cout << "Results saved to output.csv" << std::endl;
    } else {
        std::cerr << "Unable to open output.csv for writing" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;

    return EXIT_SUCCESS;
}

