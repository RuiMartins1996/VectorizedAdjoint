#ifndef CSE_HPP
#define CSE_HPP

#include <boost/numeric/odeint/integrate/check_adapter.hpp>
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

// Boost
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;

// Obtain a baseline solution for comparison
template <class System, class Stepper, class State, class Time>
std::vector<double>  cse_const(System system, Stepper stepper, const State u0, const State alphas,
           const Time ti, const Time tf,const Time dt) {
    // std::cout << "Computing baseline solution" << std::endl;
    int Nin = u0.size();
    int Npar = alphas.size();

    typedef typename odeint::unwrap_reference<Stepper>::type stepper_type;

    std::vector<double> x0aug(Nin + Nin * Npar, 0.0);
    for (int i = 0; i < Nin; i++) x0aug[i] = u0[i];

    integrate_const(stepper, system, x0aug,
        ti, tf, dt);

    std::vector<double> out(Nin * Npar);
    for (int i = 0; i < out.size(); i++) out[i] = x0aug[i + Nin];

    return out;
}

// Obtain a baseline solution for comparison
template <class System, class Stepper, class State, class Time>
std::vector<double>  cse(System system, Stepper stepper, const State u0, const State alphas,
           const Time ti, const Time tf,const Time dt, const double AbsTol,
           const double RelTol) {
    // std::cout << "Computing baseline solution" << std::endl;
    int Nin = u0.size();
    int Npar = alphas.size();

    typedef typename odeint::unwrap_reference<Stepper>::type stepper_type;

    std::vector<double> x0aug(Nin + Nin * Npar, 0.0);
    for (int i = 0; i < Nin; i++) x0aug[i] = u0[i];

    integrate_adaptive(
        odeint::make_controlled<stepper_type>(AbsTol, RelTol), system, x0aug,
        ti, tf, dt);

    std::vector<double> out(Nin * Npar);
    for (int i = 0; i < out.size(); i++) out[i] = x0aug[i + Nin];

    return out;
}

// Obtain a baseline solution for comparison
template <class System, class Stepper, class State, class Time,
          typename Observer>
std::vector<double> cse(System system, Stepper stepper, const State u0, const State alphas,
           const Time ti, const Time tf,const Time dt, const double AbsTol,
           const double RelTol, Observer observer) {
    // std::cout << "Computing baseline solution" << std::endl;
    int Nin = u0.size();
    int Npar = alphas.size();

    typedef typename odeint::unwrap_reference<Stepper>::type stepper_type;

    std::vector<double> x0aug(Nin + Nin * Npar, 0.0);
    for (int i = 0; i < Nin; i++) x0aug[i] = u0[i];

    integrate_adaptive(
        odeint::make_controlled<stepper_type>(AbsTol, RelTol), system, x0aug,
        ti, tf, dt);

    std::vector<double> out(Nin * Npar);
    for (int i = 0; i < out.size(); i++) out[i] = x0aug[i + Nin];

    observer(out, Nin);
    return out;
}

#endif