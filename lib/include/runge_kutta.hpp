#ifndef RUNGE_KUTTA_HPP_INCLUDED
#define RUNGE_KUTTA_HPP_INCLUDED

#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/throw_exception.hpp>

#include "detail/runge_kutta.hpp"

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

// With observer
template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta(
    Stepper stepper, System system, State &start_state,
    const State &alphas, Time start_time, const Time end_time,
    Time dt, void *driver, Observer observer)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return detail_runge_kutta::runge_kutta(
        stepper, system, start_state, alphas, start_time, end_time, dt, driver,
        observer, stepper_category());
};

// Without observer
template <class Stepper, class System, class State, class Time>
size_t runge_kutta(
    Stepper stepper, System system, State &start_state,
    const State &alphas, Time start_time, const Time end_time,
    Time dt, void *driver)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return detail_runge_kutta::runge_kutta(
        stepper, system, start_state, alphas, start_time, end_time, dt, driver,
        boost::numeric::odeint::null_observer(), stepper_category());
};

#endif
