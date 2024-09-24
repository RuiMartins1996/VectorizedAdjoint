#ifndef BACKPROPAGATION_HPP_INCLUDED
#define BACKPROPAGATION_HPP_INCLUDED

#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/throw_exception.hpp>
#include <boost/type_traits/is_same.hpp>
#include <stdexcept>

#include "backpropagation.hpp"

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

namespace backpropagation {

template <class Stepper, class System, class State, class Time, class Observer>
size_t backpropagation(Stepper stepper, System system, State &start_state,
                       State &parameters, Time start_time, Time end_time,
                       Time dt, Time AbsTol, Time RelTol, Observer observer) {
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;
    return detail::backpropagation(stepper, system, start_state, parameters,
                                   start_time, end_time, dt, AbsTol, RelTol,
                                   observer, stepper_category());
}

template <class Stepper, class System, class State, class Time>
size_t backpropagation(Stepper stepper, System system, State &start_state,
                       State &parameters, Time start_time, Time end_time,
                       Time dt, Time AbsTol, Time RelTol) {
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;
    return detail::backpropagation(stepper, system, start_state, parameters,
                                   start_time, end_time, dt, AbsTol, RelTol,
                                   boost::numeric::odeint::null_observer(),
                                   stepper_category());
}





}  // namespace backpropagation

#endif