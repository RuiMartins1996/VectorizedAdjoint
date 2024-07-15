#ifndef CASA_HPP
#define CASA_HPP

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

#include "AadDataNew.hpp"

namespace ublas = boost::numeric::ublas;

template <class SystemActive, class System>
class AugmentedSystem
{
private:
    System system;

    std::unique_ptr<AadDataNew> p_aad_data;

    std::vector<double> alphas;

    int Nin;
    int Npar;

public:
    AugmentedSystem(SystemActive _systemActive, System _system, const int _Nin, const int _Npar, std::vector<double> _alphas) : system(_system),
                                                                                                                                p_aad_data(std::make_unique<AadDataNew>(_Nin, _Nin, _Npar, _systemActive)),
                                                                                                                                alphas(_alphas),
                                                                                                                                Nin(_Nin),
                                                                                                                                Npar(_Npar){};

    // Move constructor by r-value
    AugmentedSystem(AugmentedSystem &&other) : system(other.system),
                                               alphas(other.alphas),
                                               Nin(other.Nin),
                                               Npar(other.Npar)
    {
        p_aad_data = std::move(other.p_aad_data);
    };

    void AddToIntegrand(const std::vector<double> &x, double t, double dt, ublas::matrix<double> &jacobian)
    {
        std::vector<double> wtimesJx(Nin);
        std::vector<double> wtimesJalpha(Npar);

        std::vector<double> u(Nin);
        for (int i = 0; i < Nin; i++)
            u[i] = x[i];

        std::vector<double> a(Nin);
        for (int i = 0; i < Nin; i++)
        {
            std::vector<double> a(Nin);
            for (int j = 0; j < Nin; j++)
                a[j] = x[Nin + i * Nin + j];

            p_aad_data->VectorTimesJacobians(a, u, alphas, t, wtimesJx, wtimesJalpha);

            for (int k = 0; k < Npar; k++)
            {
                jacobian(i, k) -= wtimesJalpha[k] * dt;
            }
        }
    }

    void rhs(const std::vector<double> &x, std::vector<double> &dxdt, double t)
    {
        // rhs returns -F(u,alpha,t)!

        // First entries of dxdt are the original system solved backwards (Nin entries)
        std::vector<double> u(Nin);
        for (int i = 0; i < Nin; i++)
            u[i] = x[i];
        std::vector<double> dudt(Nin);

        system(u, dudt, t);
        for (int i = 0; i < Nin; i++)
            dxdt[i] = dudt[i];

        // Other entries of dxdt are the adjoint equations (Nin*Nin entries)

        std::vector<double> wtimesJx(Nin);
        std::vector<double> wtimesJalpha(Npar);

        for (int i = 0; i < Nin; i++)
        {
            std::vector<double> a(Nin);
            for (int j = 0; j < Nin; j++)
                a[j] = x[Nin + i * Nin + j];

            p_aad_data->VectorTimesJacobians(a, u, alphas, t, wtimesJx, wtimesJalpha);

            for (int j = 0; j < Nin; j++)
                dxdt[Nin + i * Nin + j] = -wtimesJx[j];
        }
    }
};

// Obtain a baseline solution for comparison
template <class SystemActive, class System, class Stepper, class State, class Time>
ublas::matrix<double> casa_const(SystemActive systemActive, System system, Stepper stepper, State u0, const State alphas,
                                 const Time ti, const Time tf, Time dt)
{
    using namespace boost::numeric::odeint;
    // std::cout << "Computing baseline solution" << std::endl;
    int Nin = u0.size();
    int Npar = alphas.size();

    typedef typename odeint::unwrap_reference<Stepper>::type stepper_type;

    // Compute forward solution
    integrate_const(stepper, system, u0, ti, tf, dt);

    State uend = u0;

    State x(Nin + Nin * Nin, 0.0);
    for (int i = 0; i < Nin; i++)
        x[i] = uend[i];
    for (int i = 0; i < Nin; i++)
        x[Nin + i * Nin + i] = 1.0;

    State dxdt(Nin + Nin * Nin);

    // We wish to solve the system backwards in time, from tf to ti.
    // In boost odeint we can integrate backwards in time by setting dt to be negative and switching ti and tf

    AugmentedSystem<SystemActive, System> augmentedsystem(systemActive, system, Nin, Npar, alphas);

    auto systemfun = std::bind(
        &AugmentedSystem<SystemActive, System>::rhs,
        &augmentedsystem,
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

    using namespace boost::numeric::odeint;

    typename odeint::unwrap_reference<Stepper>::type &st = stepper;

    // Compute duidalpha_k integral with trapezoidal rule
    // int_a^b f(x)dx = sum_{k=1}^N (f(x_{k-1})+f(x{k}))/2*dx_{k}
    ublas::matrix<double> jacobian(Nin, Npar);

    size_t count = 0;

    Time end_time = tf;
    Time step = -dt;
    
    //!NEED TO CHECK TRAPEZOIDAL RULE AGAIN!
    int totalNumOfSteps = static_cast<int>(end_time / dt);
    augmentedsystem.AddToIntegrand(x, end_time, step / 2, jacobian);
    while (detail::less_with_sign(end_time, ti,
                                  step))
    {
        std::cout << count << "of " << totalNumOfSteps << "\n";
        if (boost::numeric::odeint::detail::less_with_sign(
                ti, static_cast<Time>(end_time + step), step))
        {
            step = ti - end_time;
        }

        st.do_step(systemfun, x, end_time, dt);
        end_time += step;

        augmentedsystem.AddToIntegrand(x, end_time, step, jacobian);
        ++count;
    }

    // Missing adding the edge terms f(x_N) and f(x_0) to trapezoidal rule
    augmentedsystem.AddToIntegrand(x, end_time, step / 2, jacobian);

    return jacobian;
}

// Obtain a baseline solution for comparison
template <class SystemActive, class System, class Stepper, class State, class Time>
ublas::matrix<double> casa(SystemActive systemActive, System system, Stepper stepper, State u0, const State alphas,
                           const Time ti, const Time tf, const Time dt, const double AbsTol,
                           const double RelTol)
{
    // std::cout << "Computing baseline solution" << std::endl;
    int Nin = u0.size();
    int Npar = alphas.size();

    typedef typename odeint::unwrap_reference<Stepper>::type stepper_type;

    // Compute forward solution
    integrate_adaptive(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        system, u0, ti, tf, dt);

    State x(Nin + Nin * Nin, 0.0);
    for (int i = 0; i < Nin; i++)
        x[i] = u0[i];
    for (int i = 0; i < Nin; i++)
        x[Nin + i * Nin + i] = 1.0;

    State xlast = x;
    double timelast = tf;

    // We wish to solve the system backwards in time, from tf to ti.
    // In boost odeint we can integrate backwards in time by setting dt to be negative and switching ti and tf

    AugmentedSystem<SystemActive, System> augmentedsystem(systemActive, system, Nin, Npar, alphas);

    auto systemfun = std::bind(
        &AugmentedSystem<SystemActive, System>::rhs,
        &augmentedsystem,
        std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

    auto stepper2 =
        boost::numeric::odeint::make_controlled(AbsTol, RelTol, stepper_type());

    boost::numeric::odeint::controlled_step_result res;
    boost::numeric::odeint::failed_step_checker
        fail_checker; // to throw a runtime_error if step size adjustment fails
    size_t count = 0;

    Time end_time = tf;
    Time step = -dt;

    // Compute duidalpha_k integral with trapezoidal rule
    // int_a^b f(x)dx = sum_{k=1}^N (f(x_{k-1})+f(x{k}))/2*dx_{k}
    ublas::matrix<double> jacobian(Nin, Npar);
    for(int i = 0;i<Nin;i++){
        for(int k = 0;k<Npar;k++) jacobian(i,k) = 0.0;
    }

    // Do first step
    do
    {
        res = stepper2.try_step(systemfun, x, end_time, step);
        fail_checker(); // check number of failed steps
    } while (res == boost::numeric::odeint::fail);
    augmentedsystem.AddToIntegrand(xlast, timelast, step / 2, jacobian);

    // Do remaining steps
    while (boost::numeric::odeint::detail::less_with_sign(end_time, ti,
                                                          step))
    {
        if (boost::numeric::odeint::detail::less_with_sign(
                ti, static_cast<Time>(end_time + step), step))
        {
            step = ti - end_time;
        }

        boost::numeric::odeint::controlled_step_result res;
        xlast = x;
        timelast = end_time;
        augmentedsystem.AddToIntegrand(xlast, timelast, step / 2, jacobian); // step and end_time change from first call to second call of this funtion
        // We are adding contributions from the time step at the right and the time step at the left
        // Need to copy xlast = x tho.
        do
        {
            res = stepper2.try_step(systemfun, x, end_time, step);
            fail_checker(); // check number of failed steps
        } while (res == boost::numeric::odeint::fail);

        augmentedsystem.AddToIntegrand(xlast, timelast, step / 2, jacobian);

        fail_checker.reset(); // if we reach here, the step was successful ->
                              // reset fail checker

        ++count;
    }
    // Missing contribution from edge terms f(x_0) to trapezoidal rule
    augmentedsystem.AddToIntegrand(x, end_time, step / 2, jacobian);

    return jacobian;
}
#endif