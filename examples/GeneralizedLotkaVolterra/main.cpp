#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

// Boost includes
#include <boost/numeric/odeint.hpp>

// Filesystem and chrono includes
#include <filesystem>
#include <chrono>

// MyLib includes
#include "../../lib/include/backpropagation_detail.hpp"
#include "../../lib/include/runge_kutta.hpp"

using namespace boost::numeric::odeint;

namespace fs = std::filesystem;

fs::path getExecutablePath() {
    // Initialize with a safe default relative to the executable
    fs::path result = fs::current_path();

    // Attempt to retrieve the actual executable path from the OS
    std::error_code ec;
    result = fs::canonical(fs::weakly_canonical(result, ec), ec);

    // Check if canonical path retrieval succeeded
    if (ec) {
        std::cerr << "Error retrieving executable path: " << ec.message() << std::endl;
    }

    return result;
}

std::vector<double> readAlphasFromFile(int N) 
{
    std::vector<double> alphas;
    

    // Get the path to the directory containing the executable
    fs::path exePath = getExecutablePath();

    // Construct the full path to the target folder and file
    fs::path targetFolder = exePath.parent_path() / "data" / ("N" + std::to_string(N));
    fs::path filePath = targetFolder / "alphasfile_cpp.csv";

    //std::string targetFolder = directory + "/N" + std::to_string(N);
    //std::string filePath = targetFolder + "/alphasfile_cpp.csv";

    // Check if the target folder exists
    if (!fs::exists(targetFolder) || !fs::is_directory(targetFolder)) {
        std::cerr << "Directory " << targetFolder << " does not exist or is not a directory." << std::endl;
        return alphas;
    }

    // Check if the file exists
    if (!fs::exists(filePath)) {
        std::cerr << "File " << filePath << " does not exist." << std::endl;
        return alphas;
    }

    // Open the file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file " << filePath << std::endl;
        return alphas;
    }

    // Read the contents of the file into the vector
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            try {
                alphas.push_back(std::stod(value));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid value found in file: " << value << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Value out of range in file: " << value << std::endl;
            }
        }
    }

    file.close();
    return alphas;
}

const double tol = 1e-5;

// Define the type of stepper to use to solve the system of ODEs
typedef std::vector<double> state_type;
typedef odeint::runge_kutta_cash_karp54<state_type> stepper_type;
typedef runge_kutta_cash_karp54<std::vector<idouble>> i_error_stepper_type;

// System function for the Generalized Lotka Volterra system of ODEs
template <typename T>
class LotkaVolterra
{
public:
    LotkaVolterra() = default;

    // Overload operator() to define the rhs of the system of ODEs
    void operator()(const std::vector<T> &x, std::vector<T> &dxdt,
                    const std::vector<T> &parameters, T t)
    {
        int N = x.size();

        for (int i = 0; i < x.size(); i++)
        {
            T sum = 0.0;
            for (int j = 0; j < x.size(); j++)
            {
                int id = N * (i + 1) + j;
                sum += parameters[id] * x[j];
            }
            dxdt[i] = x[i] * (parameters[i] + sum);
        }
    }
};


int main(int argc, char* argv[]) {

    // Check if the correct number of arguments is provided
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <tolerance> <N>" << std::endl;
        return 1;
    }

    // Parse tolerance from command line arguments
    double tol;
    int N;

    try {
        tol = std::stod(argv[1]);  // Convert tolerance to double
        N = std::stoi(argv[2]);    // Convert N to int
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "Argument out of range: " << e.what() << std::endl;
        return 1;
    }

    // Print the parsed values
    std::cout << "Tolerance: " << tol << std::endl;
    std::cout << "N: " << N << std::endl;

    LotkaVolterra<double> lotkaVolterra;

    std::vector<double> x0(N);
    for(int i = 0; i < N; i++) x0[i] = 0.1;
    std::vector<double> dxdt(N);
    
    // Number of parameters
    int Npar = N*N+N;


    std::vector<double> alphas = readAlphasFromFile(N);

    /*
    std::cout << "Initial conditions: x0 = [";
    for(int i = 0; i < N; i++) std::cout << x0[i] << ", ";
    std::cout << "]" << std::endl;
    */

    i_error_stepper_type rk;

    double ti = 0.0;
    double tf = 10.0;
    double dt = 1e-3;
    double AbsTol = tol;
    double RelTol = tol;

    // Create a stepper object
    stepper_type stepper;

    void *driver_handle = new Driver(stepper, LotkaVolterra<idouble>(), N, N, Npar, x0);

    auto trki = std::chrono::high_resolution_clock::now();
    runge_kutta(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        LotkaVolterra<double>(),
        x0, alphas, ti, tf, dt,driver_handle);
    
    auto trkf = std::chrono::high_resolution_clock::now();

    auto start = std::chrono::time_point_cast<std::chrono::microseconds>(trki).time_since_epoch().count();
    auto end = std::chrono::time_point_cast<std::chrono::microseconds>(trkf).time_since_epoch().count();

    auto times_microseconds_rk = end - start;
    
    auto tadji = std::chrono::high_resolution_clock::now();
    ublas::matrix<double> J = 
        backpropagation::detail::compute_jacobian(driver_handle, alphas);
    auto tadjf = std::chrono::high_resolution_clock::now();

    start = std::chrono::time_point_cast<std::chrono::microseconds>(tadji).time_since_epoch().count();
    end = std::chrono::time_point_cast<std::chrono::microseconds>(tadjf).time_since_epoch().count();

    auto times_microseconds_adj = end - start;


    std::vector<double> muAdjoints(N * Npar);
    for (int i = 0; i < N; i++) {
        for (int k = 0;k<Npar;k++) muAdjoints[i*Npar+k] = J(i,k);
    }

    std::cout << "Time forward integration: " << times_microseconds_rk << " microseconds" << std::endl;
    std::cout << "Time adjoint integration: " << times_microseconds_adj << " microseconds" << std::endl;

    return EXIT_SUCCESS;
}

