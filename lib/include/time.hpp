#ifndef TIME_INCLUDED
#define TIME_INCLUDED

#include <boost/filesystem.hpp>
#include <chrono>
#include <ctime>
#include <fstream>
#include <vector>

#include "aadcintegrate.hpp"
#include "backpropagation_detail.hpp"
#include "cse.hpp"
#include "divided_differences.hpp"
#include "runge_kutta.hpp"

template <typename Value>
double compute_mean(std::vector<Value> &v)
{
    Value mean = 0.0;
    for (int i = 0; i < v.size(); i++)
        mean += v[i];

    return static_cast<double>(mean / v.size());
};

template <typename Value>
std::pair<double, double> compute_mean_and_std_dev(std::vector<Value> &v)
{

    Value sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = static_cast<double>(sum / v.size());

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),[mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    
    if(v.size()>1){
        return std::make_pair(mean,std::sqrt(sq_sum / (v.size()-1)));
    }else{
        return std::make_pair(mean,0.0);
    }
};

class Timer
{
public:
    Timer(){};

    void Start()
    {
        m_startTimePoint = std::chrono::high_resolution_clock::now();
    };

    int Stop()
    {
        auto endTimePoint = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_startTimePoint).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimePoint).time_since_epoch().count();

        auto duration = end - start;

        return duration;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startTimePoint;
};

template <typename Function, typename... Args>
std::pair<double, double> get_execution_time(int Trep, double RelTol, Function func, Args &&...args)
{
    std::vector<int> times_microseconds(Trep);

    // Warm up iterations
    for (int i = 0; i < static_cast<int>(Trep); i++)
    {
        func(std::forward<Args>(args)...);
    }

    for (int i = 0; i < Trep; i++)
    {
        auto ti = std::chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        auto tf = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(ti).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(tf).time_since_epoch().count();

        times_microseconds[i] = end - start;
    }

    std::pair<double, double> result = compute_mean_and_std_dev(times_microseconds);

    while (static_cast<double>(result.second / result.first) > RelTol)
    {
        auto ti = std::chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        auto tf = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(ti).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(tf).time_since_epoch().count();

        times_microseconds.push_back(end - start);
        result = compute_mean_and_std_dev(times_microseconds);

        std::cout << result.first << "," << result.second << "\n";
    }

    return result;
};

template <typename Function, typename... Args>
std::pair<double, double> get_execution_time(int Trep, Function func, Args &&...args)
{
    std::vector<int> times_microseconds(Trep);

    // Warm up iterations
    for (int i = 0; i < static_cast<int>(Trep / 3) + 1; i++)
    {
        func(std::forward<Args>(args)...);
    }

    for (int i = 0; i < Trep; i++)
    {
        auto ti = std::chrono::high_resolution_clock::now();
        func(std::forward<Args>(args)...);
        auto tf = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(ti).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(tf).time_since_epoch().count();

        times_microseconds[i] = end - start;
    }

    return compute_mean_and_std_dev(times_microseconds);
};

struct TimesBP
{
public:
    std::vector<double> mean_times;
    std::vector<double> std_devs;
    TimesBP(std::vector<double> _mean_times, std::vector<double> _std_devs)
        : mean_times(_mean_times), std_devs(_std_devs){};
};

struct Times
{
public:
    double mean_time;
    double std_dev;
    Times(double _mean_time, double _std_dev)
        : mean_time(_mean_time), std_dev(_std_dev){};
};

double compute_std_deviation(std::vector<double> &v)
{
    int N = v.size();

    double v_mean = compute_mean(v);

    double sum = 0.0;

    for (int i = 0; i < N; i++)
    {
        sum += (v[i] - v_mean) * (v[i] - v_mean);
    }

    if (N > 1)
    {
        return std::sqrt(sum / static_cast<double>(N - 1));
    }
    else
    {
        return std::sqrt(sum / static_cast<double>(N));
    }
};

template <class System, class Stepper, class State, class Value>
TimesBP time_aadc(System system, Stepper stepper, State &start_state,
                  const State &parameters, Value ti, Value tf, Value dt, int Trep)
{
    double start;
    double end;

    std::vector<double> times_rec(Trep);
    std::vector<double> times_fwd(Trep);
    std::vector<double> times_bp(Trep);

    State uf(start_state.size());
    std::vector<double> wtimesJx(start_state.size());
    std::vector<double> wtimesJalpha(parameters.size());

    for (int i = 0; i < Trep; i++)
    {
        start = clock();
        AadDataForRK test(start_state.size(), parameters.size(), system, ti, tf,
                          dt);
        end = clock();

        times_rec[i] = (double(end - start) / CLOCKS_PER_SEC) * 1000;

        start = clock();
        test.Run(start_state, uf, parameters);
        end = clock();
        times_fwd[i] = (double(end - start) / CLOCKS_PER_SEC) * 1000;

        start = clock();
        std::vector<double> w(start_state.size(), 0.0);
        for (int i = 0; i < start_state.size(); i++)
        {
            w[i] = 1.0;
            test.VectorTimesJacobians(w, start_state, parameters, wtimesJx,
                                      wtimesJalpha);
            w[i] = 0.0;
        };
        end = clock();
        times_bp[i] = (double(end - start) / CLOCKS_PER_SEC) * 1000;
    };

    double mean_rec = compute_mean(times_rec);
    double std_dev_rec = compute_std_deviation(times_rec, mean_rec);

    double mean_fwd = compute_mean(times_fwd);
    double std_dev_fwd = compute_std_deviation(times_fwd, mean_fwd);

    double mean_bp = compute_mean(times_bp);
    double std_dev_bp = compute_std_deviation(times_bp, mean_bp);

    return TimesBP(std::vector<double>({mean_rec, mean_fwd, mean_bp}),
                   std::vector<double>({std_dev_rec, std_dev_fwd, std_dev_bp}));
}

template <class System, class Stepper, class State, class Value>
Times time_dd(System system, Stepper stepper, State &start_state,
              const State &parameters, Value ti, Value tf, Value dt, double AbsTol,
              double RelTol, int Trep)
{
    double start;
    double end;

    double bump = 1e-5;

    std::vector<double> times(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        State u0 = start_state;

        start = clock();

        divided_differences(stepper, system, u0, parameters, bump, ti, tf, dt,
                            AbsTol, RelTol);
        end = clock();
        double duration = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        if (i > Trep - 1)
        {
            times[i - Trep] = duration;
        }
    }

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);

    return Times(mean, std_dev);
};

template <class System, class Stepper, class State, class Value>
Times time_cse(System system, Stepper stepper, State &start_state,
               const State &parameters, Value ti, Value tf, Value dt, double AbsTol,
               double RelTol, int Trep)
{
    double start;
    double end;

    std::vector<double> times(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        State u0 = start_state;
        start = clock();
        cse(system, stepper, u0, parameters, ti, tf, dt, AbsTol, RelTol);
        end = clock();
        double duration = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        if (i > Trep - 1)
        {
            times[i - Trep] = duration;
        }
    }

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);

    return Times(mean, std_dev);
};

template <class System, class Stepper, class State>
Times time_record(System system, Stepper stepper, State &start_state,
                  const State &parameters, int Trep)
{
    double start;
    double end;

    std::vector<double> times(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        start = clock();
        void *driver_handle;
        create_driver_handle(stepper, system, start_state.size(),
                             start_state.size(), parameters.size(), start_state,
                             driver_handle);
        end = clock();
        double duration = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        if (i > Trep - 1)
        {
            times[i - Trep] = duration;
        }
        // delete_driver_handle(driver_handle);
    }

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);

    return Times(mean, std_dev);
};

template <class iSystem, class System, class Stepper, class State, class Value>
Times time_rec_over_exec(iSystem isystem, System system, Stepper stepper,
                         State &start_state, const State &parameters, Value ti,
                         Value tf, Value dt, double AbsTol, double RelTol,
                         int Trep)
{
    double start;
    double end;

    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;
    // typedef typename boost::numeric::odeint::unwrap_reference<System>::type
    //     system_type;

    // system_type sys(parameters);
    std::vector<double> times_rec(Trep);
    std::vector<double> times_exec(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        State u0 = start_state;
        start = clock();
        void *driver_handle;
        create_driver_handle(stepper, isystem, start_state.size(),
                             start_state.size(), parameters.size(), u0,
                             driver_handle);
        end = clock();

        times_rec[i] = (double(end - start) / CLOCKS_PER_SEC) * 1000;

        start = clock();
        runge_kutta(boost::numeric::odeint::make_controlled<stepper_type>(
                        AbsTol, RelTol),
                    system, u0, parameters, ti, tf, dt, driver_handle);

        backpropagation::detail::backpropagation(driver_handle, parameters);
        end = clock();

        if (i > Trep - 1)
        {
            times_exec[i - Trep] = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        };
        // delete_driver_handle(driver_handle);
    }

    std::vector<double> times(Trep);
    for (int i = 0; i < Trep; i++)
        times[i] = times_rec[i] / times_exec[i];

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);
    return Times(mean, std_dev);
};

template <class System, class Stepper, class State, class Value>
Times time_bp_no_driver(System system, Stepper stepper,
                        State &start_state, const State &parameters, Value ti, Value tf,
                        Value dt, double AbsTol, double RelTol, int Trep, void *driver_handle)
{
    double start;
    double end;

    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;
    // typedef typename boost::numeric::odeint::unwrap_reference<System>::type
    //     system_type;

    std::vector<double> times(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        State u0 = start_state;
        start = clock();

        runge_kutta(boost::numeric::odeint::make_controlled<stepper_type>(
                        AbsTol, RelTol),
                    system, u0, parameters, ti, tf, dt, driver_handle);

        backpropagation::detail::backpropagation(driver_handle, parameters);
        end = clock();
        double duration = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        if (i > Trep - 1)
        {
            times[i - Trep] = duration;
        }

        // delete_driver_handle(driver_handle);
    }

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);
    return Times(mean, std_dev);
};

template <class iSystem, class System, class Stepper, class State, class Value>
Times time_bp(iSystem isystem, System system, Stepper stepper,
              State &start_state, const State &parameters, Value ti, Value tf,
              Value dt, double AbsTol, double RelTol, int Trep)
{
    double start;
    double end;

    typedef typename boost::numeric::odeint::unwrap_reference<Stepper>::type
        stepper_type;
    // typedef typename boost::numeric::odeint::unwrap_reference<System>::type
    //     system_type;

    std::vector<double> times(Trep);

    for (int i = 0; i < Trep + Trep; i++)
    {
        State u0 = start_state;
        start = clock();
        void *driver_handle;
        create_driver_handle(stepper, isystem, start_state.size(),
                             start_state.size(), parameters.size(), u0,
                             driver_handle);
        runge_kutta(boost::numeric::odeint::make_controlled<stepper_type>(
                        AbsTol, RelTol),
                    system, u0, parameters, ti, tf, dt, driver_handle);

        backpropagation::detail::backpropagation(driver_handle, parameters);
        end = clock();
        double duration = (double(end - start) / CLOCKS_PER_SEC) * 1000;
        if (i > Trep - 1)
        {
            times[i - Trep] = duration;
        }

        // delete_driver_handle(driver_handle);
    }

    double mean = compute_mean(times);
    double std_dev = compute_std_deviation(times, mean);
    return Times(mean, std_dev);
};

#endif
