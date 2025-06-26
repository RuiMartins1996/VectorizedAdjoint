#ifndef DETAIL_RUNGE_KUTTA_HPP_INCLUDED
#define DETAIL_RUNGE_KUTTA_HPP_INCLUDED

#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>

#include "../Driver.hpp"

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

namespace vectorizedadjoint
{
namespace detail_runge_kutta
{

template <class System, class State, class Time>
class SystemFunctor
{
  private:
    System system;
    State parameters;

  public:
    SystemFunctor(System _system, State _parameters) : system(_system), parameters(_parameters){};

    void operator()(const State &u, State &dudt, const Time &t)
    {
        this->system(u, dudt, parameters, t);
    };
};

// for non-controlled steppers
template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta(
    Stepper stepper, System system, State &start_state,
    const State &alphas, Time start_time,
    const Time end_time, Time dt, Driver &driver,
    Observer observer, boost::numeric::odeint::stepper_tag)
{
    typename odeint::unwrap_reference<Observer>::type &obs = observer;
    typename odeint::unwrap_reference<Stepper>::type &st = stepper;

    using namespace boost::numeric::odeint;

    auto system_runge_kutta_step =
        SystemFunctor<System, State, Time>(system, alphas);

    Time time = start_time;
    int step = 0;
    // cast time+dt explicitely in case of expression templates (e.g. multiprecision)
    while (boost::numeric::odeint::detail::less_eq_with_sign(static_cast<Time>(time + dt), end_time, dt)) {
        obs(start_state, time);
        driver.PushBackTime(time);
        driver.PushBackState(start_state);

        st.do_step(system_runge_kutta_step, start_state, time, dt);
        // direct computation of the time avoids error propagation happening when using time += dt
        // we need clumsy type analysis to get boost units working here
        ++step;
        time = start_time + static_cast<typename unit_value_type<Time>::type>(step) * dt;
    }

    obs(start_state, time);
    driver.PushBackTime(time);
    driver.PushBackState(start_state);

    return step;
}

// For controlled steppers
template <class Stepper, class System, class State, class Time, class Observer>
size_t runge_kutta(
    Stepper stepper, System system, State &start_state,
    const State &alphas, Time start_time,
    const Time end_time, Time dt, Driver &driver,
    Observer observer, odeint::controlled_stepper_tag)
{
    typename odeint::unwrap_reference<Observer>::type &obs = observer;
    typename boost::numeric::odeint::unwrap_reference<Stepper>::type &st = stepper;

    boost::numeric::odeint::failed_step_checker
        fail_checker; // to throw a runtime_error if step size adjustment fails
    size_t count = 0;

    auto system_runge_kutta_step =
        SystemFunctor<System, State, Time>(system, alphas);

    driver.p_states->Clear();
    while (boost::numeric::odeint::detail::less_with_sign(start_time, end_time, dt)) {
        obs(start_state, start_time);
        driver.PushBackTime(start_time);
        driver.PushBackState(start_state);

        if (boost::numeric::odeint::detail::less_with_sign(
                end_time, static_cast<Time>(start_time + dt), dt)) {
            dt = end_time - start_time;
        }

        boost::numeric::odeint::controlled_step_result res;
        do {
            res = st.try_step(system_runge_kutta_step, start_state, start_time, dt);
            fail_checker(); // check number of failed steps
        } while (res == boost::numeric::odeint::fail);
        fail_checker.reset(); // if we reach here, the step was successful ->
                              // reset fail checker
        ++count;
    }

    obs(start_state, start_time);
    driver.PushBackTime(start_time);
    driver.PushBackState(start_state);

    return count;
}

} // namespace detail_runge_kutta
} // namespace vectorizedadjoint
#endif