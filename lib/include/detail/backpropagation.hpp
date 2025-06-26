#ifndef DETAIL_BACKPROPAGATION_HPP_INCLUDED
#define DETAIL_BACKPROPAGATION_HPP_INCLUDED

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
#include <boost/type_traits/is_same.hpp>

// #include "OdeDriver.hpp"
#include "../Driver.hpp"

namespace vectorizedadjoint
{
namespace detail
{

template <class State, class Time>
ublas::matrix<double> compute_intermediate_states(Driver &driver, int Nin, const State &alphas, const Time time, const Time dt, int n)
{
    int s = driver.p_butcher->GetStages();

    ublas::matrix<double> K_vectors(Nin, s + 1);
    State dxdt(Nin);

    State u(Nin);

    driver.GetState(u, n);

    for (int m = 0; m < s; m++) {
        // Get intermediate state x^{m,n}
        State xmn(Nin);

        xmn = u;

        for (int j = 0; j < m; j++)
            for (int i = 0; i < Nin; i++)
                xmn[i] += dt * driver.p_butcher->a(m, j) * K_vectors(i, j);

        State dxdt(Nin);

        driver.Rhs(xmn, dxdt, alphas, time);

        for (int i = 0; i < Nin; i++)
            K_vectors(i, m) = dxdt[i];
    };

    // Store dudt in last row of K_vectors
    State dudt(Nin);
    for (int i = 0; i < Nin; i++) {
        dudt[i] = 0.0;
        for (int j = 0; j < s; j++)
            dudt[i] += driver.p_butcher->b(j) * K_vectors(i, j);
        K_vectors(i, s) = u[i] + dt * dudt[i];
    };

    return K_vectors;
};

template <class State, class Time>
State get_intermediate_state(Driver &driver, const State &u, const ublas::matrix<double> &K_vectors_n, const Time dt, const int m, const int n)
{
    int Nin = K_vectors_n.size1();

    //(k=0)
    State u_mn = u;

    //(k=1...m)
    for (int k = 1; k < m; k++) {
        for (int i = 0; i < Nin; i++)
            u_mn[i] += dt * driver.p_butcher->a(m - 1, k - 1) * K_vectors_n(i, k - 1);
    }

    return u_mn;
}

template <class State>
void back_prop_step(Driver &driver, aadc::mmVector<mmType> &wbarend, aadc::mmVector<mmType> &alphabar, const State &alphas, int n)
{
    int Nin = wbarend.size();
    int Npar = alphas.size();

    double time = driver.GetTime(n);
    double dt = driver.GetDt(n);

    int s = driver.p_butcher->GetStages();

    // Initialize w_bar_n
    std::vector<aadc::mmVector<mmType>> w_bar_n(s + 2);
    for (int m = 0; m < w_bar_n.size(); m++) {
        aadc::mmVector<mmType> temp(Nin);
        for (int i = 0; i < Nin; i++)
            temp[i] = aadc::mmZero<mmType>();
        w_bar_n[m] = temp;
    }

    // Link adjoints at different iterations
    w_bar_n[s + 1] = wbarend;

    // Step s+1
    for (int i = 0; i < Nin; i++) {
        w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], w_bar_n[s + 1][i]);
        for (int m = 1; m < s + 1; m++) {
            w_bar_n[m][i] = aadc::mmAdd(w_bar_n[m][i], driver.p_butcher->b(m - 1) * dt * w_bar_n[s + 1][i]);
        }
    }

    std::vector<double> wk(Nin);
    aadc::mmVector<mmType> w_bar_m_n(Nin);
    aadc::mmVector<mmType> wJx(Nin);
    aadc::mmVector<mmType> wJalpha(Npar);

    // Get u(t^n)
    State u(Nin);
    driver.GetState(u, n);
    // And recompute this step
    ublas::matrix<double> K_vectors_n = compute_intermediate_states(driver, Nin, alphas, time, dt, n);

    // Steps s to 1
    for (int m = s; m > 0; m--) {
        double time_m_n = time + driver.p_butcher->c(m) * dt;

        //* Get the mth intermediate x^mn state and F derivatives
        wk = get_intermediate_state(driver, u, K_vectors_n, dt, m, n);

        //* Get values of (m)th adjoint on temporary matrix
        for (int i = 0; i < Nin; i++)
            w_bar_m_n[i] = w_bar_n[m][i];

        //* Do vector times Jacobian multiplication (vA) for state and
        driver.p_aad_data->VectorTimesJacobians(w_bar_m_n, wk, alphas, time_m_n, wJx, wJalpha);
        //* Update adjoints of order 0...m-1
        for (int i = 0; i < Nin; i++) {
            // Update m=0 adjoint
            w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], wJx[i]);
            // Update kth<mth adjoint
            for (int k = 1; k < m; k++) {
                w_bar_n[k][i] = aadc::mmAdd(w_bar_n[k][i], wJx[i] * driver.p_butcher->a(m - 1, k - 1) * dt);
            }
        };

        //* Update alphabar
        if (m > 0) {
            for (int k = 0; k < Npar; k++)
                alphabar[k] = aadc::mmAdd(alphabar[k], wJalpha[k]);
        }
    };

    //* Update wbarend for next iteration
    for (int i = 0; i < Nin; i++)
        wbarend[i] = w_bar_n[0][i];
}

template <class State>
void back_prop_step(Driver &driver, State &wbarend, State &alphabar, const State &alphas, int n)
{
    int Nin = wbarend.size();
    int Npar = alphas.size();

    double time = driver.GetTime(n);
    double dt = driver.GetDt(n);

    int s = driver.p_butcher->GetStages();
    // Initialize w_bar_n
    ublas::matrix<double> w_bar_n(Nin, s + 2);
    for (int i = 0; i < Nin; i++)
        for (int j = 0; j < w_bar_n.size2() - 1; j++)
            w_bar_n(i, j) = 0.0;
    for (int i = 0; i < Nin; i++)
        w_bar_n(i, s + 1) = wbarend[i];

    // Step s+1
    for (int i = 0; i < Nin; i++) {
        w_bar_n(i, 0) += w_bar_n(i, s + 1);
        for (int m = 1; m < s + 1; m++) {
            w_bar_n(i, m) += driver.p_butcher->b(m - 1) * dt * w_bar_n(i, s + 1);
        }
    }
    std::vector<double> wk(Nin);
    std::vector<double> w_bar_m_vec(Nin);
    std::vector<double> wJx(Nin);
    std::vector<double> wJalpha(Npar);

    // Get u(t^n)
    State u(Nin);
    driver.GetState(u, n);
    // And recompute this step
    ublas::matrix<double> K_vectors_n = compute_intermediate_states(driver, Nin, alphas, time, dt, n);
    // Steps s to 1
    for (int m = s; m > 0; m--) {
        double time_m_n = time + driver.p_butcher->c(m) * dt;

        //* Get the mth intermediate x^mn state and F derivatives
        wk = get_intermediate_state(driver, u, K_vectors_n, dt, m, n);

        //* Get values of (m)th adjoint on temporary matrix
        for (int i = 0; i < Nin; i++)
            w_bar_m_vec[i] = w_bar_n(i, m);

        //* Do vector times Jacobian multiplication (vA) for state and parameter Jacobians
        driver.p_aad_data->VectorTimesJacobians(w_bar_m_vec, wk, alphas, time_m_n, wJx, wJalpha);

        //* Update adjoints of order 0...m-1
        for (int i = 0; i < Nin; i++) {
            // Update m=0 adjoint
            w_bar_n(i, 0) += wJx[i];
            // Update kth<mth adjoint
            for (int k = 1; k < m; k++) {
                w_bar_n(i, k) += wJx[i] * driver.p_butcher->a(m - 1, k - 1) * dt;
            }
        };

        //* Update alphabar
        if (m > 0) {
            for (int k = 0; k < Npar; k++)
                alphabar[k] += wJalpha[k];
        }
    };

    //* Update wbarend for next iteration
    for (int i = 0; i < Nin; i++)
        wbarend[i] = w_bar_n(i, 0);
}

template <class State, class Observer>
aadc::mmVector<mmType> back_prop(Driver &driver, aadc::mmVector<mmType> &wbarend, const State &alphas, Observer observer)
{
    int num_time_steps = driver.GetT() - 1;
    int Npar = alphas.size();

    // Initliaze alphabar;
    aadc::mmVector<mmType> alphabar(Npar);
    for (int k = 0; k < alphabar.size(); k++)
        alphabar[k] = aadc::mmZero<mmType>();

    double time = driver.GetTime(num_time_steps);
    observer(wbarend, time);

    //(iteration n=T-1 to n=1)
    for (int n = num_time_steps - 1; n > -1; n--) {
        back_prop_step(driver, wbarend, alphabar, alphas, n);

        time = driver.GetTime(n);
        observer(wbarend, time);
    }

    return alphabar;
};

template <class State, class Observer>
State back_prop(Driver &driver, State &wbarend, const State &alphas, Observer observer)
{
    int num_time_steps = driver.GetT() - 1;
    int Npar = alphas.size();

    // Initliaze alphabar;
    State alphabar(Npar);
    for (int k = 0; k < alphabar.size(); k++)
        alphabar[k] = 0.0;

    double time = driver.GetTime(num_time_steps);
    observer(wbarend, time);
    //(iteration n=T-1 to n=1)
    for (int n = num_time_steps - 1; n > -1; n--) {
        back_prop_step(driver, wbarend, alphabar, alphas, n);

        time = driver.GetTime(n);
        observer(wbarend, time);
    }

    return alphabar;
};

template <class State>
void adjointSolve(Driver &driver, const State &parameters)
{
    // This function is going to change the entries of the vectors driver.p_lambda and driver.p_mu!
    int Nout = driver.GetNout();
    int Npar = driver.GetNpar();
    int Nin = driver.GetNin();

    // Number of doubles that fit into an mmType. We can "process" a batch at a time.
    int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));

    int numOfBatches = Nout / sizeOfBatch;
    int remainderEntries = Nout % sizeOfBatch;

    aadc::mmVector<mmType> wbarend(Nin);
    aadc::mmVector<mmType> alphabar(Npar);

    // Handle reverse pass for first numOfBatches*sizeOfBatch objetive functions
    // Each batch of sizeOfBatch objective functions is "processed" in parallel
    for (int n = 0; n < numOfBatches; n++) {

        // Load seed with the values of lambda
        for (int j = 0; j < sizeOfBatch; j++) { // Iterate over the vectors in a single batch
            for (int i = 0; i < Nin; i++) {     // Fill the seed with the values of lambda
                wbarend[i][j] = (*driver.p_lambda)[n * sizeOfBatch + j][i];
            }
        }

        alphabar = back_prop(driver, wbarend, parameters, boost::numeric::odeint::null_observer());

        // Update lambda
        for (int j = 0; j < sizeOfBatch; j++) { // Iterate over the vectors in a single batch
            for (int i = 0; i < Nin; i++) {     // Fill the seed with the values of lambda
                (*driver.p_lambda)[n * sizeOfBatch + j][i] = wbarend[i][j];
            }
        }

        for (int j = 0; j < sizeOfBatch; j++) { // Iterate over the vectors in a single batch
            for (int k = 0; k < Npar; k++) {
                (*driver.p_mu)[n * sizeOfBatch + j][k] += alphabar[k][j];
            }
        }
    }

    // Handle reverse pass for last Nout % sizeOfBatch objetive functions

    // Load seed with the values of lambda
    for (int j = 0; j < remainderEntries; j++) { // Iterate over the vectors in a single batch
        for (int i = 0; i < Nin; i++) {          // Fill the seed with the values of lambda
            wbarend[i][j] = (*driver.p_lambda)[(numOfBatches)*sizeOfBatch + j][i];
        }
    }

    alphabar = back_prop(driver, wbarend, parameters, boost::numeric::odeint::null_observer());

    // Update lambda
    // TODO: There must be a better way to do this than to copy the values to wbarend from lambda and the copy the result wbarend to lambda
    for (int j = 0; j < remainderEntries; j++) { // Iterate over the vectors in a single batch
        for (int i = 0; i < Nin; i++) {          // Fill the seed with the values of lambda
            (*driver.p_lambda)[(numOfBatches)*sizeOfBatch + j][i] = wbarend[i][j];
        }
    }

    for (int j = 0; j < remainderEntries; j++) { // Iterate over the vectors in a single batch
        for (int k = 0; k < Npar; k++) {
            (*driver.p_mu)[(numOfBatches)*sizeOfBatch + j][k] += alphabar[k][j];
        }
    }
};

//! LEGACY:

// This function handles the computation of the sensitivity matrix of the ODE system WITHOUT leveraging SIMD vectorization
template <class State>
auto computeSensitivityMatrixNoSIMD(Driver &driver, const State &parameters)
{
    int Nin = driver.GetNin();
    int Npar = driver.GetNpar();

    std::vector<std::vector<double>> adjoints(Nin, std::vector<double>(Npar, 0.0));

    for (int i = 0; i < Nin; i++) {
        State wbarend(Nin, 0.0);
        wbarend[i] = 1.0;

        State alphabar = back_prop(driver, wbarend, parameters, boost::numeric::odeint::null_observer());
        for (int j = 0; j < alphabar.size(); j++)
            adjoints[i][j] = alphabar[j];
    }

    return adjoints;
};

// This function handles the computation of the sensitivity matrix of the ODE system leveraging SIMD vectorization
template <class State>
auto computeSensitivityMatrix(Driver &driver, const State &parameters)
{
    int Nin = driver.GetNin();
    int Npar = driver.GetNpar();

    std::vector<std::vector<double>> jacobian(Nin, std::vector<double>(Npar, 0.0));

    aadc::mmVector<mmType> wbarend(Nin);
    aadc::mmVector<mmType> alphabar(Npar);

    // Lambda
    auto resetSeed = [](aadc::mmVector<mmType> &wbarend) {
        for (int i = 0; i < wbarend.size(); i++) {
            wbarend[i] = aadc::mmZero<mmType>();
        }
    };

    int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));

    int numOfBatches = Nin / sizeOfBatch;
    int remainderEntries = Nin % sizeOfBatch;

    // Handle entries of wbarend below maximum multiple of sizeOfBatch
    for (int n = 0; n < numOfBatches; n++) {
        // std::cout << "Batch " << n <<" of "<< numOfBatches <<" \n";
        //  Reset Seed
        resetSeed(wbarend);

        // Load seed
        for (int j = 0; j < sizeOfBatch; j++)
            wbarend[n * sizeOfBatch + j][j] = 1.0;

        alphabar = back_prop(driver, wbarend, parameters, boost::numeric::odeint::null_observer());

        for (int j = 0; j < sizeOfBatch; j++) {
            for (int k = 0; k < parameters.size(); k++) {
                jacobian[n * sizeOfBatch + j][k] = alphabar[k][j];
            }
        }
    }

    // Handle remaining entries

    // Reset Seed
    resetSeed(wbarend);

    // Load seed
    for (int j = 0; j < remainderEntries; j++)
        wbarend[(numOfBatches)*sizeOfBatch + j][j] = 1.0;

    alphabar = back_prop(driver, wbarend, parameters, boost::numeric::odeint::null_observer());

    for (int j = 0; j < remainderEntries; j++) {
        for (int k = 0; k < parameters.size(); k++) {
            jacobian[(numOfBatches)*sizeOfBatch + j][k] = alphabar[k][j];
        }
    }

    return jacobian;
};

} // namespace detail
} // namespace vectorizedadjoint

#endif
