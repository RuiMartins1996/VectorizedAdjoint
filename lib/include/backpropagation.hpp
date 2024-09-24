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

template <class State>
std::vector<std::vector<double>> backpropagation(void *driver, const State &parameters)
{
    return detail::backpropagation(driver, parameters);
};

template <class State>
ublas::matrix<double> compute_jacobian(void *driver, const State &parameters)
{
    return detail::compute_jacobian(driver, parameters);
};

// With observer
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

} // end namespace backpropagation

#endif
