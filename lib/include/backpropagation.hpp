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
#include "OdeDriverNew.hpp"
#include "utilities.hpp"

// For parallelization
#include <future>
#include "ThreadPool.hpp"
#include <chrono>

static std::mutex mtx;

namespace backpropagation
{
    namespace detail
    {

        void write_to_jacobian_from_thread(
            ublas::matrix<double> &jacobian,
            aadc::mmVector<mmType> &alphabar,
            int n)
        {
            int Nin = jacobian.size1();
            int Npar = jacobian.size2();

            int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));
            int numOfBatches = Nin / sizeOfBatch;
            int remainderEntries = Nin % sizeOfBatch;

            // std::cout << "Batch number:" << n << "\n";
            if (n < numOfBatches)
            {
                for (int j = 0; j < sizeOfBatch; j++)
                {
                    for (int k = 0; k < Npar; k++)
                    {
                        jacobian(n * sizeOfBatch + j, k) = alphabar[k][j];
                        // std::cout<< "Writing in:("<<n * sizeOfBatch + j<<","<<k<<")"<< std::endl;
                    }
                }
            }
            else
            {
                for (int j = 0; j < remainderEntries; j++)
                {
                    for (int k = 0; k < Npar; k++)
                    {
                        jacobian((numOfBatches)*sizeOfBatch + j, k) = alphabar[k][j];
                        // std::cout<< "Writing in:("<<(numOfBatches)*sizeOfBatch + j<<","<<k<<")"<< std::endl;
                    }
                }
            }

            // log(jacobian);
            // }
        };

        template <class State, class Time>
        ublas::matrix<double> compute_intermediate_states(Driver *p_driver, int Nin,
                                                          const State &alphas,
                                                          const Time time,
                                                          const Time dt, int n)
        {
            int s = p_driver->p_butcher->GetStages();

            ublas::matrix<double> K_vectors(Nin, s + 1);
            State dxdt(Nin);

            State u(Nin);

            p_driver->GetState(u, n);

            for (int m = 0; m < s; m++)
            {
                // Get intermediate state x^{m,n}
                State xmn(Nin);

                xmn = u;

                for (int j = 0; j < m; j++)
                {
                    for (int i = 0; i < Nin; i++)
                        xmn[i] += dt * p_driver->p_butcher->a(m, j) * K_vectors(i, j);
                }

                State dxdt(Nin);

                p_driver->Rhs(xmn, dxdt, alphas, time);

                for (int i = 0; i < Nin; i++)
                {
                    K_vectors(i, m) = dxdt[i];
                }
            };

            // Store dudt in last row of K_vectors
            State dudt(Nin);
            for (int i = 0; i < Nin; i++)
            {
                dudt[i] = 0.0;
                for (int j = 0; j < s; j++)
                    dudt[i] += p_driver->p_butcher->b(j) * K_vectors(i, j);
                K_vectors(i, s) = u[i] + dt * dudt[i];
            };

            return K_vectors;
        };

        template <class State, class Time>
        State get_intermediate_state(Driver *p_driver, const State &u,
                                     const ublas::matrix<double> &K_vectors_n,
                                     const Time dt, const int m, const int n)
        {
            int Nin = K_vectors_n.size1();

            //(k=0)
            State u_mn = u;

            //(k=1...m)
            for (int k = 1; k < m; k++)
            {
                for (int i = 0; i < Nin; i++)
                    u_mn[i] += dt * p_driver->p_butcher->a(m - 1, k - 1) *
                               K_vectors_n(i, k - 1);
            }

            return u_mn;
        }

        template <class State>
        void back_prop_step_parallel(Driver *p_driver, aadc::mmVector<mmType> &wbarend,
                                     aadc::mmVector<mmType> &alphabar, const State &alphas,
                                     int n, const int numOfBatch, ublas::matrix<double> &jacobian)
        {

            int Nin = wbarend.size();
            int Npar = alphas.size();

            int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));
            int numOfBatches = Nin / sizeOfBatch;
            int remainderEntries = Nin % sizeOfBatch;

            double time = p_driver->GetTime(n);
            double dt = p_driver->GetDt(n);

            int s = p_driver->p_butcher->GetStages();

            // Initialize w_bar_n
            std::vector<aadc::mmVector<mmType>> w_bar_n(s + 2);
            for (int m = 0; m < w_bar_n.size(); m++)
            {
                aadc::mmVector<mmType> temp(Nin);
                for (int i = 0; i < Nin; i++)
                    temp[i] = aadc::mmZero<mmType>();
                w_bar_n[m] = temp;
            }

            // Link adjoints at different iterations
            w_bar_n[s + 1] = wbarend;

            // Step s+1
            for (int i = 0; i < Nin; i++)
            {
                w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], w_bar_n[s + 1][i]);
                for (int m = 1; m < s + 1; m++)
                {
                    w_bar_n[m][i] =
                        aadc::mmAdd(w_bar_n[m][i], p_driver->p_butcher->b(m - 1) * dt *
                                                       w_bar_n[s + 1][i]);
                }
            }

            std::vector<double> wk(Nin);
            aadc::mmVector<mmType> w_bar_m_n(Nin);
            aadc::mmVector<mmType> wJx(Nin);
            aadc::mmVector<mmType> wJalpha(Npar);

            // Get u(t^n)
            State u(Nin);
            p_driver->GetState(u, n);
            // And recompute this step
            mtx.lock(); // prevent data race
            ublas::matrix<double> K_vectors_n =
                compute_intermediate_states(p_driver, Nin, alphas, time, dt, n);
            mtx.unlock();

            // Steps s to 1
            for (int m = s; m > 0; m--)
            {
                double time_m_n = time + p_driver->p_butcher->c(m) * dt;

                //* Get the mth intermediate x^mn state and F derivatives
                wk = get_intermediate_state(p_driver, u, K_vectors_n, dt, m, n);

                //* Get values of (m)th adjoint on temporary matrix
                for (int i = 0; i < Nin; i++)
                    w_bar_m_n[i] = w_bar_n[m][i];

                //* Do vector times Jacobian multiplication (vA) for state and
                mtx.lock(); // prevent data race
                p_driver->p_aad_data->VectorTimesJacobians(w_bar_m_n, wk, alphas,
                                                           time_m_n, wJx, wJalpha);
                mtx.unlock();

                //* Update adjoints of order 0...m-1
                for (int i = 0; i < Nin; i++)
                {
                    // Update m=0 adjoint
                    w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], wJx[i]);
                    // Update kth<mth adjoint
                    for (int k = 1; k < m; k++)
                    {
                        w_bar_n[k][i] = aadc::mmAdd(
                            w_bar_n[k][i],
                            wJx[i] * p_driver->p_butcher->a(m - 1, k - 1) * dt);
                    }
                };

                //* Update alphabar
                if (m > 0)
                {
                    if (n < numOfBatches)
                    {
                        for (int j = 0; j < sizeOfBatch; j++)
                        {
                            for (int k = 0; k < Npar; k++)
                            {
                                jacobian(n * sizeOfBatch + j, k) += wJalpha[k][j];
                                //std::cout<< "Writing in:("<<n * sizeOfBatch + j<<","<<k<<")"<< std::endl;
                            }
                        }
                    }
                    else
                    {
                        for (int j = 0; j < remainderEntries; j++)
                        {
                            for (int k = 0; k < Npar; k++)
                            {
                                jacobian((numOfBatches)*sizeOfBatch + j, k) += wJalpha[k][j];
                                //std::cout<< "Writing in:("<<(numOfBatches)*sizeOfBatch + j<<","<<k<<")"<< std::endl;
                            }
                        }
                    }
                    /*
                    for (int k = 0; k < Npar; k++)
                        alphabar[k] = aadc::mmAdd(alphabar[k], wJalpha[k]);
                    */
                }
            };

            //* Update wbarend for next iteration
            for (int i = 0; i < Nin; i++)
                wbarend[i] = w_bar_n[0][i];
        }

        template <class State>
        void back_prop_step(Driver *p_driver, aadc::mmVector<mmType> &wbarend,
                            aadc::mmVector<mmType> &alphabar, const State &alphas,
                            int n)
        {

            int Nin = wbarend.size();
            int Npar = alphas.size();

            double time = p_driver->GetTime(n);
            double dt = p_driver->GetDt(n);

            int s = p_driver->p_butcher->GetStages();

            // Initialize w_bar_n
            std::vector<aadc::mmVector<mmType>> w_bar_n(s + 2);
            for (int m = 0; m < w_bar_n.size(); m++)
            {
                aadc::mmVector<mmType> temp(Nin);
                for (int i = 0; i < Nin; i++)
                    temp[i] = aadc::mmZero<mmType>();
                w_bar_n[m] = temp;
            }

            // Link adjoints at different iterations
            w_bar_n[s + 1] = wbarend;

            // Step s+1
            for (int i = 0; i < Nin; i++)
            {
                w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], w_bar_n[s + 1][i]);
                for (int m = 1; m < s + 1; m++)
                {
                    w_bar_n[m][i] =
                        aadc::mmAdd(w_bar_n[m][i], p_driver->p_butcher->b(m - 1) * dt *
                                                       w_bar_n[s + 1][i]);
                }
            }

            std::vector<double> wk(Nin);
            aadc::mmVector<mmType> w_bar_m_n(Nin);
            aadc::mmVector<mmType> wJx(Nin);
            aadc::mmVector<mmType> wJalpha(Npar);

            // Get u(t^n)
            State u(Nin);
            p_driver->GetState(u, n);
            // And recompute this step
            mtx.lock(); // prevent data race
            ublas::matrix<double> K_vectors_n =
                compute_intermediate_states(p_driver, Nin, alphas, time, dt, n);
            mtx.unlock();

            // Steps s to 1
            for (int m = s; m > 0; m--)
            {
                double time_m_n = time + p_driver->p_butcher->c(m) * dt;

                //* Get the mth intermediate x^mn state and F derivatives
                wk = get_intermediate_state(p_driver, u, K_vectors_n, dt, m, n);

                //* Get values of (m)th adjoint on temporary matrix
                for (int i = 0; i < Nin; i++)
                    w_bar_m_n[i] = w_bar_n[m][i];

                //* Do vector times Jacobian multiplication (vA) for state and
                mtx.lock(); // prevent data race
                p_driver->p_aad_data->VectorTimesJacobians(w_bar_m_n, wk, alphas,
                                                           time_m_n, wJx, wJalpha);
                mtx.unlock();

                //* Update adjoints of order 0...m-1
                for (int i = 0; i < Nin; i++)
                {
                    // Update m=0 adjoint
                    w_bar_n[0][i] = aadc::mmAdd(w_bar_n[0][i], wJx[i]);
                    // Update kth<mth adjoint
                    for (int k = 1; k < m; k++)
                    {
                        w_bar_n[k][i] = aadc::mmAdd(
                            w_bar_n[k][i],
                            wJx[i] * p_driver->p_butcher->a(m - 1, k - 1) * dt);
                    }
                };

                //* Update alphabar
                if (m > 0)
                {
                    for (int k = 0; k < Npar; k++)
                        alphabar[k] = aadc::mmAdd(alphabar[k], wJalpha[k]);
                }
            };

            //* Update wbarend for next iteration
            for (int i = 0; i < Nin; i++)
                wbarend[i] = w_bar_n[0][i];
        }

        template <class State>
        void back_prop_step(Driver *p_driver, State &wbarend, State &alphabar,
                            const State &alphas, int n)
        {
            int Nin = wbarend.size();
            int Npar = alphas.size();

            double time = p_driver->GetTime(n);
            double dt = p_driver->GetDt(n);

            int s = p_driver->p_butcher->GetStages();
            // Initialize w_bar_n
            ublas::matrix<double> w_bar_n(Nin, s + 2);
            for (int i = 0; i < Nin; i++)
                for (int j = 0; j < w_bar_n.size2() - 1; j++)
                    w_bar_n(i, j) = 0.0;
            for (int i = 0; i < Nin; i++)
                w_bar_n(i, s + 1) = wbarend[i];

            // Step s+1
            for (int i = 0; i < Nin; i++)
            {
                w_bar_n(i, 0) += w_bar_n(i, s + 1);
                for (int m = 1; m < s + 1; m++)
                {
                    w_bar_n(i, m) +=
                        p_driver->p_butcher->b(m - 1) * dt * w_bar_n(i, s + 1);
                }
            }
            std::vector<double> wk(Nin);
            std::vector<double> w_bar_m_vec(Nin);
            std::vector<double> wJx(Nin);
            std::vector<double> wJalpha(Npar);

            // Get u(t^n)
            State u(Nin);
            p_driver->GetState(u, n);
            // And recompute this step
            ublas::matrix<double> K_vectors_n =
                compute_intermediate_states(p_driver, Nin, alphas, time, dt, n);
            // Steps s to 1
            for (int m = s; m > 0; m--)
            {
                double time_m_n = time + p_driver->p_butcher->c(m) * dt;

                //* Get the mth intermediate x^mn state and F derivatives
                wk = get_intermediate_state(p_driver, u, K_vectors_n, dt, m, n);

                //* Get values of (m)th adjoint on temporary matrix
                for (int i = 0; i < Nin; i++)
                    w_bar_m_vec[i] = w_bar_n(i, m);

                //* Do vector times Jacobian multiplication (vA) for state and
                // parameter Jacobians

                // vector_times_jacobian(p_driver,w_bar_m_vec, wk, alphas,
                //                                           time, wJx, wJalpha);

                p_driver->p_aad_data->VectorTimesJacobians(w_bar_m_vec, wk, alphas,
                                                           time_m_n, wJx, wJalpha);

                //(lldb) breakpoint set --name breakpoint --condition 'time < 0.01' ;

                //* Update adjoints of order 0...m-1
                for (int i = 0; i < Nin; i++)
                {
                    // Update m=0 adjoint
                    w_bar_n(i, 0) += wJx[i];
                    // Update kth<mth adjoint
                    for (int k = 1; k < m; k++)
                    {
                        w_bar_n(i, k) +=
                            wJx[i] * p_driver->p_butcher->a(m - 1, k - 1) * dt;
                    }
                };

                //* Update alphabar
                if (m > 0)
                {
                    for (int k = 0; k < Npar; k++)
                        alphabar[k] += wJalpha[k];
                }
            };

            //* Update wbarend for next iteration
            for (int i = 0; i < Nin; i++)
                wbarend[i] = w_bar_n(i, 0);
        }

        template <class State>
        void back_prop_parallel(Driver *p_driver,
                                aadc::mmVector<mmType> &wbarend,
                                const State &alphas, ublas::matrix<double> &jacobian, int numOfBatch)
        {
            int num_time_steps = p_driver->GetT() - 1;
            int Npar = alphas.size();

            // Initliaze alphabar;
            aadc::mmVector<mmType> alphabar(Npar);
            for (int k = 0; k < alphabar.size(); k++)
                alphabar[k] = aadc::mmZero<mmType>();

            //(iteration n=T-1 to n=1)
            for (int n = num_time_steps - 1; n > -1; n--)
            {
                back_prop_step(p_driver, wbarend, alphabar, alphas, n);
            }

            write_to_jacobian_from_thread(jacobian, alphabar, numOfBatch);

            // return alphabar;
        };

        template <class State, class Observer>
        aadc::mmVector<mmType> back_prop(Driver *p_driver,
                                         aadc::mmVector<mmType> &wbarend,
                                         const State &alphas, Observer observer)
        {
            int num_time_steps = p_driver->GetT() - 1;
            int Npar = alphas.size();

            // Initliaze alphabar;
            aadc::mmVector<mmType> alphabar(Npar);
            for (int k = 0; k < alphabar.size(); k++)
                alphabar[k] = aadc::mmZero<mmType>();

            double time = p_driver->GetTime(num_time_steps);
            observer(wbarend, time);

            //(iteration n=T-1 to n=1)
            for (int n = num_time_steps - 1; n > -1; n--)
            {
                back_prop_step(p_driver, wbarend, alphabar, alphas, n);

                time = p_driver->GetTime(n);
                observer(wbarend, time);
            }

            return alphabar;
        };

        template <class State, class Observer>
        State back_prop(Driver *p_driver, State &wbarend, const State &alphas,
                        Observer observer)
        {
            int num_time_steps = p_driver->GetT() - 1;
            int Npar = alphas.size();

            // Initliaze alphabar;
            State alphabar(Npar);
            for (int k = 0; k < alphabar.size(); k++)
                alphabar[k] = 0.0;

            double time = p_driver->GetTime(num_time_steps);
            observer(wbarend, time);
            //(iteration n=T-1 to n=1)
            for (int n = num_time_steps - 1; n > -1; n--)
            {
                back_prop_step(p_driver, wbarend, alphabar, alphas, n);

                time = p_driver->GetTime(n);
                observer(wbarend, time);
            }

            return alphabar;
        };

        template <class State>
        State backpropagation(void *driver, State &wbarend, const State &parameters)
        {
            Driver *p_driver = static_cast<Driver *>(driver);
            return back_prop(p_driver, wbarend, parameters,
                             boost::numeric::odeint::null_observer());
        };

        template <class State, class Observer>
        State backpropagation(void *driver, State &wbarend, const State &parameters,
                              Observer observer)
        {
            Driver *p_driver = static_cast<Driver *>(driver);
            return back_prop(p_driver, wbarend, parameters, observer);
        };

        //!  USING THIS NOW
        template <class State, class Observer>
        size_t backpropagation(void *driver, const State &parameters,
                               Observer observer)
        {
            Driver *p_driver = static_cast<Driver *>(driver);

            int Nin = p_driver->GetNin();
            int Npar = p_driver->GetNpar();

            //! Just for saving correctly
            std::vector<std::vector<double>> alphasmat;
            //! Just for saving correctly
            for (int i = 0; i < Nin; i++)
            {
                std::vector<double> temp(Npar);
                alphasmat.push_back(temp);
            }

            for (int i = 0; i < Nin; i++)
            {
                State wbarend(Nin, 0.0);
                wbarend[i] = 1.0;

                State alphabar = back_prop(p_driver, wbarend, parameters,
                                           boost::numeric::odeint::null_observer());
                //! Just for saving correctly
                for (int j = 0; j < alphabar.size(); j++)
                    alphasmat[i][j] = alphabar[j];
            }

            //! Just for saving correctly
            std::vector<double> alphafinal(Nin * Npar);
            //! Just for saving correctly
            convert(alphafinal, alphasmat);

            observer(alphafinal, Nin);

            return 1;
        };

        //! USING THIS NOW
        template <class State>
        std::vector<double> backpropagation(void *driver, const State &parameters)
        {
            Driver *p_driver = static_cast<Driver *>(driver);

            int Nin = p_driver->GetNin();
            int Npar = p_driver->GetNpar();

            //! Just for saving correctly
            std::vector<std::vector<double>> alphasmat;
            //! Just for saving correctly
            for (int i = 0; i < Nin; i++)
            {
                std::vector<double> temp(Npar);
                alphasmat.push_back(temp);
            }

            for (int i = 0; i < Nin; i++)
            {
                State wbarend(Nin, 0.0);
                wbarend[i] = 1.0;

                State alphabar = back_prop(p_driver, wbarend, parameters,
                                           boost::numeric::odeint::null_observer());
                //! Just for saving correctly
                for (int j = 0; j < alphabar.size(); j++)
                    alphasmat[i][j] = alphabar[j];
            }

            //! Just for saving correctly
            std::vector<double> alphafinal(Nin * Npar);
            //! Just for saving correctly
            convert(alphafinal, alphasmat);

            return alphafinal;
        };

        template <class State>
        ublas::matrix<double> compute_jacobian(void *driver, const State &parameters)
        {
            Driver *p_driver = static_cast<Driver *>(driver);

            int Nin = p_driver->GetNin();
            int Npar = p_driver->GetNpar();

            ublas::matrix<double> jacobian(Nin, Npar);

            aadc::mmVector<mmType> wbarend(Nin);
            aadc::mmVector<mmType> alphabar(Npar);

            // Lambda
            auto resetSeed = [](aadc::mmVector<mmType> &wbarend)
            {
                for (int i = 0; i < wbarend.size(); i++)
                {
                    wbarend[i] = aadc::mmZero<mmType>();
                }
            };

            int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));

            int numOfBatches = Nin / sizeOfBatch;
            int remainderEntries = Nin % sizeOfBatch;

            // Handle entries of wbarend below maximum multiple of sizeOfBatch
            for (int n = 0; n < numOfBatches; n++)
            {   
                //std::cout << "Batch " << n <<" of "<< numOfBatches <<" \n";
                // Reset Seed
                resetSeed(wbarend);

                // Load seed
                for (int j = 0; j < sizeOfBatch; j++)
                    wbarend[n * sizeOfBatch + j][j] = 1.0;

                alphabar = back_prop(p_driver, wbarend, parameters,
                                     boost::numeric::odeint::null_observer());

                for (int j = 0; j < sizeOfBatch; j++)
                {
                    for (int k = 0; k < parameters.size(); k++)
                    {
                        jacobian(n * sizeOfBatch + j, k) = alphabar[k][j];
                    }
                }
            }

            // Handle remaining entries

            // Reset Seed
            resetSeed(wbarend);

            // Load seed
            for (int j = 0; j < remainderEntries; j++)
                wbarend[(numOfBatches)*sizeOfBatch + j][j] = 1.0;

            alphabar = back_prop(p_driver, wbarend, parameters,
                                 boost::numeric::odeint::null_observer());

            for (int j = 0; j < remainderEntries; j++)
            {
                for (int k = 0; k < parameters.size(); k++)
                {
                    jacobian((numOfBatches)*sizeOfBatch + j, k) = alphabar[k][j];
                }
            }

            /*
            //! Just for saving correctly
            std::vector<double> alphafinal(Nin * Npar);
            //! Just for saving correctly
            convert(alphafinal, jacobian);

            return alphafinal;
            */
            return jacobian;
        };

        void write_to_jacobian_from_thread(
            ublas::matrix<double> &jacobian,
            std::vector<std::future<aadc::mmVector<mmType>>> &futures,
            int n)
        {
            int Nin = jacobian.size1();
            int Npar = jacobian.size2();

            int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));
            int numOfBatches = Nin / sizeOfBatch;
            int remainderEntries = Nin % sizeOfBatch;

            // if (futures[n].valid() && futures[n].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
            //{
            aadc::mmVector<mmType> alphabar = futures[n].get();

            if (n < numOfBatches)
            {
                for (int j = 0; j < sizeOfBatch; j++)
                {
                    for (int k = 0; k < Npar; k++)
                    {
                        jacobian(n * sizeOfBatch + j, k) = alphabar[k][j];
                    }
                }
            }
            else
            {
                for (int j = 0; j < remainderEntries; j++)
                {
                    for (int k = 0; k < Npar; k++)
                    {
                        jacobian((numOfBatches)*sizeOfBatch + j, k) = alphabar[k][j];
                    }
                }
            }
            //}
        };

        // In this function we do the Jacobian computation du(tf)/dalpha and du(tf)du0
        // This function employs vectorization via operating over aadc::mmVector<mmType>
        // This function employs parallelization via using std::future
        template <class State>
        ublas::matrix<double> compute_jacobian_2(
            void *driver, const State &parameters)
        {
            Driver *p_driver = static_cast<Driver *>(driver);

            int Nin = p_driver->GetNin();
            int Npar = p_driver->GetNpar();

            ublas::matrix<double> jacobian(Nin, Npar);
            for(int i = 0;i<Nin;i++){
                for(int k = 0;k<Npar;k++) jacobian(i,k) = 0.0;
            }

            aadc::mmVector<mmType> wbarend(Nin);

            aadc::mmVector<mmType> alphabar(Npar);

            // Lambda
            auto resetSeed = [](aadc::mmVector<mmType> &u)
            {
                for (int i = 0; i < u.size(); i++)
                {
                    u[i] = aadc::mmZero<mmType>();
                }
            };

            int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));

            int numOfBatches = Nin / sizeOfBatch;
            int remainderEntries = Nin % sizeOfBatch;

            // ThreadPool thread_pool;

            // thread_pool.Start();

            // One future for each batch and one other for the remaining entries
            //!
            // std::vector<std::future<aadc::mmVector<mmType>>> futures;

            std::vector<std::future<void>> futures;
            // Handle batches
            for (int n = 0; n < numOfBatches + !(remainderEntries == 0); n++)
            {
                // Reset Seed
                resetSeed(wbarend);

                if (n > numOfBatches - 1) // Handle remaining entries
                {
                    // Load seed
                    for (int j = 0; j < remainderEntries; j++)
                        wbarend[(numOfBatches)*sizeOfBatch + j][j] = 1.0;

                    /*
                    thread_pool.QueueJob([&, wbarend,n]()
                        {
                            int _n = n;
                            aadc::mmVector<mmType> _wbarend = wbarend;
                            return back_prop_parallel(p_driver, _wbarend, parameters, std::ref(jacobian), _n);
                        });
                    */

                    futures.emplace_back(std::async(
                        std::launch::async,
                        [&, wbarend, n]()
                        {
                            int _n = n;
                            aadc::mmVector<mmType> _wbarend = wbarend;
                            return back_prop_parallel(p_driver, _wbarend, parameters, std::ref(jacobian), _n);
                            // return back_prop(p_driver, _wbarend, parameters, boost::numeric::odeint::null_observer());
                        }));
                    break;
                } // else, handle batches

                // Load seed
                for (int j = 0; j < sizeOfBatch; j++)
                    wbarend[n * sizeOfBatch + j][j] = 1.0;

                /*
                thread_pool.QueueJob([&, wbarend,n]()
                        {
                            int _n = n;
                            aadc::mmVector<mmType> _wbarend = wbarend;
                            return back_prop_parallel(p_driver, _wbarend, parameters, std::ref(jacobian), _n);
                        });
                */
                futures.emplace_back(std::async(
                    std::launch::async,
                    [&, wbarend, n]()
                    {
                        int _n = n;
                        aadc::mmVector<mmType> _wbarend = wbarend;
                        return back_prop_parallel(p_driver, _wbarend, parameters, std::ref(jacobian), _n);
                        // return back_prop(p_driver, _wbarend, parameters, boost::numeric::odeint::null_observer());
                    }));
            }

            // Wait for futures to be ready
            while (std::none_of(
                futures.begin(),
                futures.end(), [](std::future<void> &future)
                { return (future.wait_for(std::chrono::seconds(0)) == std::future_status::ready); }))
            {
            }

            // while there exist uncompleted tasks
            // iterate over uncompleted tasks
            // retrieve result to jacobian matrix

            /*
            while (std::any_of(
                futures.begin(),
                futures.end(), [](std::future<aadc::mmVector<mmType>> &future)
                { return future.valid(); })) // While there are valid futurues
            {

                // Loop through valid/ready futures and write
                for (int n = 0; n < numOfBatches + !(remainderEntries == 0); n++)
                {
                    if (futures[n].valid() && futures[n].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
                    {
                        write_to_jacobian_from_thread(jacobian, futures, n);
                    }
                }

                if (std::none_of(
                        futures.begin(),
                        futures.end(), [](std::future<aadc::mmVector<mmType>> &future)
                        { return future.valid(); }))
                {
                    break;
                }

                // Check if ALL futures are ready
                if (std::all_of(
                        futures.begin(),
                        futures.end(), [](std::future<aadc::mmVector<mmType>> &future)
                        { return future.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }))
                {
                    std::cout << "All Threads Completed" << std::endl;
                    break;
                }
            }
            */

            /*
            //! Just for saving correctly
            std::vector<double> alphafinal(Nin * Npar);
            //! Just for saving correctly
            convert(alphafinal, jacobian);

            return alphafinal;
            */
            
            return jacobian;
        };

    } // namespace detail
} // namespace backpropagation

#endif