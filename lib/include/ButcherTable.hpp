#ifndef BUTCHERTABLE_HPP_INCLUDED
#define BUTCHERTABLE_HPP_INCLUDED

// Boost
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace odeint = boost::numeric::odeint;

class ButcherTable
{

  private:
    std::unique_ptr<std::vector<double>> p_a;
    std::unique_ptr<std::vector<double>> p_b;
    std::unique_ptr<std::vector<double>> p_c;

  public:
    template <typename Stepper>
    ButcherTable(Stepper stepper)
    {
        typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category stepper_category;

        initialize(stepper, stepper_category());
    };

    template <typename Stepper>
    void initialize(Stepper stepper, odeint::explicit_controlled_stepper_tag)
    {
        std::cout << "To create Driver, we need to supply an error stepper or a regular stepper, not a controlled stepper!" << std::endl;
    }

    // For non-controlled steppers
    template <typename Stepper>
    void initialize(Stepper stepper, odeint::stepper_tag)
    {
        int stages;

        std::unique_ptr<std::vector<double>> _a;
        std::unique_ptr<std::vector<double>> _b;
        std::unique_ptr<std::vector<double>> _c;

        using State = typename odeint::unwrap_reference<Stepper>::type::state_type;

        if (typeid(stepper) == typeid(odeint::euler<State>)) {
            // std::cout << "Euler\n";

            stages = 1;

            _a = std::make_unique<std::vector<double>>(stages * stages);
            _b = std::make_unique<std::vector<double>>(stages);
            _c = std::make_unique<std::vector<double>>(stages);

            std::fill(_a->begin(), _a->end(), 0.0);
            (*_b)[0] = 1.0;
            (*_c)[0] = 0.0;

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        } else if (typeid(stepper) == typeid(odeint::runge_kutta4_classic<State>) || typeid(stepper) == typeid(odeint::runge_kutta4<State>)) {
            // std::cout << "Runge-Kutta 4\n";

            using namespace odeint;

            stages = 4;

            boost::array<double, 1> a1 = rk4_coefficients_a1<double>();
            boost::array<double, 2> a2 = rk4_coefficients_a2<double>();
            boost::array<double, 3> a3 = rk4_coefficients_a3<double>();

            // Store a1, a2, a3 in a tuple
            auto a_tuple = std::make_tuple(
                rk4_coefficients_a1<double>(),
                rk4_coefficients_a2<double>(),
                rk4_coefficients_a3<double>());

            _a = std::make_unique<std::vector<double>>(stages * stages);
            _b = std::make_unique<std::vector<double>>(stages);
            _c = std::make_unique<std::vector<double>>(stages);

            // Initialize the matrix _a to zero
            std::fill(_a->begin(), _a->end(), 0.0);

            // Lambda function to get the index of the matrix _a
            auto id = [stages](int i, int j) { return i * stages + j; };

            // Define the inner lambda that processes each array
            auto process_array = [&_a, &id](const auto &arr, int &row) {
                for (std::size_t i = 0; i < arr.size(); ++i) {
                    (*_a)[id(row, i)] = arr[i];
                }
                ++row; // Move to the next row for the next array
            };

            // Define the outer lambda, which handles the unpacking of the tuple and calling process_array for each element.
            auto process_tuple = [&_a, &id, &process_array](auto &&...args) {
                int row = 1;
                // Call process_array for each element in the args parameter pack
                (process_array(args, row), ...);
            };

            // Apply the process_tuple lambda to the a_tuple
            std::apply(process_tuple, a_tuple);

            // Fill in _b and _c vectors
            boost::array<double, 4> b = rk4_coefficients_b<double>();
            boost::array<double, 4> c = rk4_coefficients_c<double>();

            std::copy(b.begin(), b.end(), _b->begin());
            std::copy(c.begin(), c.end(), _c->begin());

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        } else {
            std::cout << "This non-controlled stepper is not supported yet!" << std::endl;
            std::terminate();
        };
    }

    // For controlled (error) steppers
    template <typename Stepper>
    void initialize(Stepper stepper, odeint::error_stepper_tag)
    {
        using namespace odeint;

        using State = typename odeint::unwrap_reference<Stepper>::type::state_type;

        int stages;

        std::unique_ptr<std::vector<double>> _a;
        std::unique_ptr<std::vector<double>> _b;
        std::unique_ptr<std::vector<double>> _c;

        if (typeid(stepper) == typeid(odeint::runge_kutta_cash_karp54<State>)) {
            // std::cout << "Runge-Kutta CashKarp(45) 4\n";

            stages = 6;

            // Store a1, a2, a3, a4, a5 in a tuple
            auto a_tuple = std::make_tuple(
                rk54_ck_coefficients_a1<double>(),
                rk54_ck_coefficients_a2<double>(),
                rk54_ck_coefficients_a3<double>(),
                rk54_ck_coefficients_a4<double>(),
                rk54_ck_coefficients_a5<double>());

            _a = std::make_unique<std::vector<double>>(stages * stages);
            _b = std::make_unique<std::vector<double>>(stages);
            _c = std::make_unique<std::vector<double>>(stages);

            // Initialize the matrix _a to zero
            std::fill(_a->begin(), _a->end(), 0.0);

            // Lambda function to get the index of the matrix _a
            auto id = [stages](int i, int j) { return i * stages + j; };

            // Define the inner lambda that processes each array
            auto process_array = [&_a, &id](const auto &arr, int &row) {
                for (std::size_t i = 0; i < arr.size(); ++i) {
                    (*_a)[id(row, i)] = arr[i];
                }
                ++row; // Move to the next row for the next array
            };

            // Define the outer lambda, which handles the unpacking of the tuple and calling process_array for each element.
            auto process_tuple = [&_a, &id, &process_array](auto &&...args) {
                int row = 1;
                // Call process_array for each element in the args parameter pack
                (process_array(args, row), ...);
            };

            // Apply the process_tuple lambda to the a_tuple
            std::apply(process_tuple, a_tuple);

            boost::array<double, 6> b = odeint::rk54_ck_coefficients_b<double>();
            boost::array<double, 6> c = odeint::rk54_ck_coefficients_c<double>();

            std::copy(b.begin(), b.end(), _b->begin());
            std::copy(c.begin(), c.end(), _c->begin());

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        } else if (typeid(stepper) == typeid(odeint::runge_kutta_fehlberg78<State>)) {
            // std::cout << "Runge-Kutta Fehlberg(78) 4\n";

            int stages = 13;

            auto a_tuple = std::make_tuple(
                rk78_coefficients_a1<double>(),
                rk78_coefficients_a2<double>(),
                rk78_coefficients_a3<double>(),
                rk78_coefficients_a4<double>(),
                rk78_coefficients_a5<double>(),
                rk78_coefficients_a6<double>(),
                rk78_coefficients_a7<double>(),
                rk78_coefficients_a8<double>(),
                rk78_coefficients_a9<double>(),
                rk78_coefficients_a10<double>(),
                rk78_coefficients_a11<double>(),
                rk78_coefficients_a12<double>());

            _a = std::make_unique<std::vector<double>>(stages * stages);
            _b = std::make_unique<std::vector<double>>(stages);
            _c = std::make_unique<std::vector<double>>(stages);

            // Initialize the matrix _a to zero
            std::fill(_a->begin(), _a->end(), 0.0);

            // Lambda function to get the index of the matrix _a
            auto id = [stages](int i, int j) { return i * stages + j; };

            // Define the inner lambda that processes each array
            auto process_array = [&_a, &id](const auto &arr, int &row) {
                for (std::size_t i = 0; i < arr.size(); ++i) {
                    (*_a)[id(row, i)] = arr[i];
                }
                ++row; // Move to the next row for the next array
            };

            // Define the outer lambda, which handles the unpacking of the tuple and calling process_array for each element.
            auto process_tuple = [&_a, &id, &process_array](auto &&...args) {
                int row = 1;
                // Call process_array for each element in the args parameter pack
                (process_array(args, row), ...);
            };

            // Apply the process_tuple lambda to the a_tuple
            std::apply(process_tuple, a_tuple);

            boost::array<double, 13> b = rk78_coefficients_b<double>();
            boost::array<double, 13> c = rk78_coefficients_c<double>();

            std::copy(b.begin(), b.end(), _b->begin());
            std::copy(c.begin(), c.end(), _c->begin());

            p_a = std::move(_a);
            p_b = std::move(_b);
            p_c = std::move(_c);
        } else {
            std::cout << "This controlled stepper is not supported yet!" << std::endl;
            std::terminate();
        };
    };

  public:
    int GetStages() { return (*p_b).size(); };
    double a(int i, int j)
    {
        int s = GetStages();
        return (*p_a)[s * i + j];
    };

    double b(int i) { return (*p_b)[i]; };

    double c(int i) { return (*p_c)[i]; }
};

#endif
