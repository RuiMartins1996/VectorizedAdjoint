#ifndef RUNGE_KUTTA_HPP_INCLUDED
#define RUNGE_KUTTA_HPP_INCLUDED

#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>
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

#include "OdeDriverNew.hpp"
#include "SystemFunctor.hpp"

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

template <class Stepper, class System, class State, class Time>
size_t runge_kutta(Stepper stepper, System system, State &start_state,
                   const State &alphas, Time start_time, const Time end_time,
                   Time dt, void *driver)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return runge_kutta_detail(stepper, system, start_state, alphas, start_time,
                              end_time, dt, driver, stepper_category());
};

template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta(Stepper stepper, System system, State &start_state,
                   const State &alphas, Time start_time, const Time end_time,
                   Time dt, void *driver, Observer observer)
{
    typedef typename odeint::unwrap_reference<Stepper>::type::stepper_category
        stepper_category;

    return runge_kutta_detail(stepper, system, start_state, alphas, start_time,
                              end_time, dt, driver, observer,
                              stepper_category());
};

template <class Stepper, class System, class State, class Time>
size_t runge_kutta_detail(Stepper stepper, System system, State &start_state,
                          const State &alphas, Time start_time,
                          const Time end_time, Time dt, void *driver,
                          odeint::controlled_stepper_tag)
{

    typename boost::numeric::odeint::unwrap_reference<Stepper>::type &st =
        stepper;

    boost::numeric::odeint::failed_step_checker
        fail_checker; // to throw a runtime_error if step size adjustment fails
    size_t count = 0;

    Driver *p_driver = static_cast<Driver *>(driver);

    auto system_runge_kutta_step =
        SystemFunctor<System, State, Time>(system, alphas);

    p_driver->p_states->Clear();
    while (boost::numeric::odeint::detail::less_with_sign(start_time, end_time,
                                                          dt))
    {
        p_driver->PushBackTime(start_time);
        p_driver->PushBackState(start_state);

        if (boost::numeric::odeint::detail::less_with_sign(
                end_time, static_cast<Time>(start_time + dt), dt))
        {
            dt = end_time - start_time;
        }

        boost::numeric::odeint::controlled_step_result res;
        do
        {
            res = st.try_step(system_runge_kutta_step, start_state, start_time,
                              dt);
            fail_checker(); // check number of failed steps
        } while (res == boost::numeric::odeint::fail);
        fail_checker.reset(); // if we reach here, the step was successful ->
                              // reset fail checker

        ++count;
    }

    p_driver->PushBackTime(start_time);
    p_driver->PushBackState(start_state);

    return count;
}

template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta_detail(Stepper stepper, System system, State &start_state,
                          const State &alphas, Time start_time,
                          const Time end_time, Time dt, void *driver,
                          Observer observer, odeint::controlled_stepper_tag)
{
    typename odeint::unwrap_reference<Observer>::type &obs = observer;
    typename boost::numeric::odeint::unwrap_reference<Stepper>::type &st =
        stepper;

    boost::numeric::odeint::failed_step_checker
        fail_checker; // to throw a runtime_error if step size adjustment fails
    size_t count = 0;

    Driver *p_driver = static_cast<Driver *>(driver);

    auto system_runge_kutta_step =
        SystemFunctor<System, State, Time>(system, alphas);

    p_driver->p_states->Clear();
    while (boost::numeric::odeint::detail::less_with_sign(start_time, end_time,
                                                          dt))
    {
        p_driver->PushBackTime(start_time);
        p_driver->PushBackState(start_state);

        if (boost::numeric::odeint::detail::less_with_sign(
                end_time, static_cast<Time>(start_time + dt), dt))
        {
            dt = end_time - start_time;
        }

        boost::numeric::odeint::controlled_step_result res;
        do
        {
            res = st.try_step(system_runge_kutta_step, start_state, start_time,
                              dt);
            fail_checker(); // check number of failed steps
        } while (res == boost::numeric::odeint::fail);
        fail_checker.reset(); // if we reach here, the step was successful ->
                              // reset fail checker

        ++count;
    }

    p_driver->PushBackTime(start_time);
    p_driver->PushBackState(start_state);
    obs(start_state, start_state.size());

    return count;
}

template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta_detail(Stepper stepper, System system, State &start_state,
                          const State &alphas, Time start_time,
                          const Time end_time, Time dt, void *driver,
                          Observer observer, boost::numeric::odeint::stepper_tag)
{

    using namespace boost::numeric::odeint;
    Driver *p_driver = static_cast<Driver *>(driver);

    typename odeint::unwrap_reference<Observer>::type &obs = observer;
    typename odeint::unwrap_reference<Stepper>::type &st = stepper;

    Time time = start_time;
    int step = 0;
    // cast time+dt explicitely in case of expression templates (e.g. multiprecision)
    while (detail::less_eq_with_sign(static_cast<Time>(time + dt), end_time, dt))
    {
        obs(start_state, time);
        p_driver->PushBackTime(time);
        p_driver->PushBackState(start_state);
        st.do_step(system, start_state, time, dt);
        // direct computation of the time avoids error propagation happening when using time += dt
        // we need clumsy type analysis to get boost units working here
        ++step;
        time = start_time + static_cast<typename unit_value_type<Time>::type>(step) * dt;
    }
    
    obs(start_state, time);
    p_driver->PushBackTime(time);
    p_driver->PushBackState(start_state);

    return step;
}

#endif
