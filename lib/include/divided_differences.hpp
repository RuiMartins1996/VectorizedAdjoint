#ifndef DIVIDED_DIFFERENCES_HPP
#define DIVIDED_DIFFERENCES_HPP

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
#include <ctime>


namespace ublas = boost::numeric::ublas;


void convertvector(std::vector<double> &alphafinal,ublas::matrix<double> &adjoints){
    int Nin = adjoints.size1();
    int Npar = adjoints.size2();

    for (int i = 0;i<Nin;i++){
        for (int k = 0;k<Npar;k++){
            int id = Npar*i+k;
            alphafinal[id] = adjoints(i,k);
        } 
    }
}




template <class Stepper, class System,class ObjectiveFunction, class State, class Time>
std::vector<double> divided_differences(Stepper stepper, System system,ObjectiveFunction objective,
                                        State &start_state, const State &parameters,
                                        double bump, Time start_time,
                                        Time end_time, Time dt, double AbsTol,
                                        double RelTol) {
    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;

    State parameters_bumped = parameters;

    ublas::matrix<double> adjoints(1,parameters.size());

    double bump_1 = 1.0 / bump;

    State u = start_state;
    // We define the parameters of the system with operator []
    system.set_parameters(parameters);

    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        system, u, start_time, end_time, dt);

    std::vector<double> res1(1);
    objective(u,res1,parameters,end_time);

    for (int k = 0; k < parameters.size(); k++) {
        State u_bumped = start_state;

        parameters_bumped[k] = parameters_bumped[k] + bump;

        system.set_parameters(parameters_bumped);

        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled<stepper_type>(AbsTol,
                                                                  RelTol),
            system, u_bumped, start_time, end_time, dt);

        std::vector<double> res2(1);
        objective(u_bumped,res2,parameters_bumped,end_time);


        for (int i = 0; i < adjoints.size1(); i++) {
            double value = (res2[i] - res1[i]) * bump_1;
            adjoints(i, k) = value;
        }

        parameters_bumped[k] = parameters[k];
    };

    std::vector<double> alphafinal(1 * parameters.size());
    for (int i = 0;i<parameters.size();i++) alphafinal[i] = adjoints(0,i);


    return alphafinal;
};


template <class Stepper, class System, class State, class Time>
std::vector<double> divided_differences_const(Stepper stepper, System system,
                                        State &start_state, const State &parameters,
                                        double bump, Time start_time,
                                        Time end_time, Time dt) {
    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;

    State parameters_bumped = parameters;

    ublas::matrix<double> adjoints(start_state.size(), parameters.size());

    double bump_1 = 1.0 / bump;

    State u = start_state;
    // We define the parameters of the system with operator []
    system.set_parameters(parameters);

    boost::numeric::odeint::integrate_const(
        stepper,system, u, start_time, end_time, dt);

    for (int k = 0; k < parameters.size(); k++) {
        State u_bumped = start_state;

        parameters_bumped[k] = parameters_bumped[k] + bump;

        system.set_parameters(parameters_bumped);

        boost::numeric::odeint::integrate_const(
            stepper,
            system, u_bumped, start_time, end_time, dt);

        for (int i = 0; i < adjoints.size1(); i++) {
            adjoints(i, k) = (u_bumped[i] - u[i]) * bump_1;
        }

        parameters_bumped[k] = parameters[k];
    };

    //! Just for saving correctly
    std::vector<double> alphafinal(start_state.size() * parameters.size());
    //! Just for saving correctly
    convertvector(alphafinal, adjoints);

    return alphafinal;
};


template <class Stepper, class System, class State, class Time>
std::vector<double> divided_differences(Stepper stepper, System system,
                                        State &start_state, const State &parameters,
                                        double bump, Time start_time,
                                        Time end_time, Time dt, double AbsTol,
                                        double RelTol) {
    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;

    State parameters_bumped = parameters;

    ublas::matrix<double> adjoints(start_state.size(), parameters.size());

    double bump_1 = 1.0 / bump;

    State u = start_state;
    // We define the parameters of the system with operator []
    system.set_parameters(parameters);

    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        system, u, start_time, end_time, dt);

    for (int k = 0; k < parameters.size(); k++) {
        State u_bumped = start_state;

        parameters_bumped[k] = parameters_bumped[k] + bump;

        system.set_parameters(parameters_bumped);

        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled<stepper_type>(AbsTol,
                                                                  RelTol),
            system, u_bumped, start_time, end_time, dt);

        for (int i = 0; i < adjoints.size1(); i++) {
            adjoints(i, k) = (u_bumped[i] - u[i]) * bump_1;
        }

        parameters_bumped[k] = parameters[k];
    };

    //! Just for saving correctly
    std::vector<double> alphafinal(start_state.size() * parameters.size());
    //! Just for saving correctly
    convertvector(alphafinal, adjoints);

    return alphafinal;
};

template <class Stepper, class System, class State, class Time, class Observer>
std::vector<double> divided_differences(Stepper stepper, System system,
                                        State &start_state, const State &parameters,
                                        double bump, Time start_time,
                                        Time end_time, Time dt, double AbsTol,
                                        double RelTol, Observer observer) {
       typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;
    State parameters_bumped = parameters;

    ublas::matrix<double> adjoints(start_state.size(), parameters.size());

    double bump_1 = 1.0 / bump;

    State u = start_state;
    // We define the parameters of the system with operator []
    system.set_parameters(parameters);

    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::make_controlled<stepper_type>(AbsTol, RelTol),
        system, u, start_time, end_time, dt);

    for (int k = 0; k < parameters.size(); k++) {
        State u_bumped = start_state;

        parameters_bumped[k] = parameters_bumped[k] + bump;

        system.set_parameters(parameters_bumped);

        boost::numeric::odeint::integrate_adaptive(
            boost::numeric::odeint::make_controlled<stepper_type>(AbsTol,
                                                                  RelTol),
            system, u_bumped, start_time, end_time, dt);

        for (int i = 0; i < adjoints.size1(); i++) {
            adjoints(i, k) = (u_bumped[i] - u[i]) * bump_1;
        }

        parameters_bumped[k] = parameters[k];
    };

    //! Just for saving correctly
    std::vector<double> alphafinal(start_state.size() * parameters.size());
    //! Just for saving correctly
   // convert(alphafinal, adjoints);

    observer(alphafinal, start_state.size());

    return alphafinal;
};

#endif