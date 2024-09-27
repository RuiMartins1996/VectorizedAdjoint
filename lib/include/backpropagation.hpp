#ifndef BACKPROPAGATION_HPP_INCLUDED
#define BACKPROPAGATION_HPP_INCLUDED

#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include "detail/backpropagation.hpp"

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

namespace backpropagation
{

// Computes the adjoints of the cost functions (the total derivatives of the cost function with respect to the initial conditions and with respect to the parameters)
// TODO: This does not use the SIMD vectorization yet
template <class State>
void adjointSolve(Driver &driver, const State &parameters)
{

    try {
        if (!(driver.p_lambda && driver.p_mu)) {
            // Either p_lambda and p_mu are null
            throw std::runtime_error("Must call setCostGradients() first!");
        }

        if (!driver.p_butcher) {
            // p_butcher is not null
            throw std::runtime_error("Must call constructDriverButcherTableau() to set Butcher Tableau!");
        }

        if (!driver.p_aad_data) {
            // aad_data is not null
            throw std::runtime_error("Must call recordDriverRHSFunction() to record the RHS with automatic differentiation!");
        }

        detail::adjointSolve(driver, parameters);

    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
    }
};

// Computes the sensitivity matrix of the ODE system WITHOUT SIMD vectorization
template <class State>
auto computeSensitivityMatrixNoSIMD(Driver &driver, const State &parameters)
{
    try {
        if (!driver.p_butcher) {
            // p_butcher is not null
            throw std::runtime_error("Must call constructDriverButcherTableau() to set Butcher Tableau!");
        }

        if (!driver.p_aad_data) {
            // aad_data is not null
            throw std::runtime_error("Must call recordDriverRHSFunction() to record the RHS with automatic differentiation!");
        }

        return detail::computeSensitivityMatrixNoSIMD(driver, parameters);
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return std::vector(0, std::vector<double>(0));
    }
};

// Computes the sensitivity matrix of the ODE system leveraging SIMD vectorization
template <class State>
auto computeSensitivityMatrix(Driver &driver, const State &parameters)
{
    try {
        if (!driver.p_butcher) {
            // p_butcher is not null
            throw std::runtime_error("Must call constructDriverButcherTableau() to set Butcher Tableau!");
        }

        if (!driver.p_aad_data) {
            // aad_data is not null
            throw std::runtime_error("Must call recordDriverRHSFunction() to record the RHS with automatic differentiation!");
        }

        return detail::computeSensitivityMatrix(driver, parameters);
    } catch (std::exception &e) {
        std::cout << e.what() << std::endl;
        return std::vector(0, std::vector<double>(0));
    }
};

// With observer
/*
template <class Stepper, class System, class State, class Time, class Observer>
size_t backpropagation(
    Stepper stepper, System system, State &start_state,
    State &parameters, Time start_time, Time end_time,
    Time dt, Time AbsTol, Time RelTol, Observer observer)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return detail::backpropagation(
        stepper, system, start_state, parameters,
        start_time, end_time, dt, AbsTol, RelTol,
        observer, stepper_category());
}

// Without observer
template <class Stepper, class System, class State, class Time>
size_t backpropagation(
    Stepper stepper, System system, State &start_state,
    State &parameters, Time start_time, Time end_time,
    Time dt, Time AbsTol, Time RelTol)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return detail::backpropagation(
        stepper, system, start_state, parameters,
        start_time, end_time, dt, AbsTol, RelTol,
        boost::numeric::odeint::null_observer(),
        stepper_category());
}
*/

} // end namespace backpropagation

#endif
