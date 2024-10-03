#ifndef ODEDRIVER_NEW_HPP // guard
#define ODEDRIVER_NEW_HPP

#include "AadData.hpp"
#include "ButcherTable.hpp"
#include "StateStorage.hpp"

// Boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace lib
{

class Driver
{
  public:
    std::unique_ptr<AadData> p_aad_data;
    std::unique_ptr<StateStorage> p_states;
    std::unique_ptr<ButcherTable> p_butcher;

    std::vector<std::vector<double>> *p_lambda;
    std::vector<std::vector<double>> *p_mu;

  private:
    int m_Nin;
    int m_Nout;
    int m_Npar;

  public:
    template <class Stepper, class System>
    Driver(Stepper stepper, System system, int Nin, int Nout, int Npar) : p_aad_data(std::make_unique<AadData>(Nin, Nout, Npar, system)),
                                                                          p_states(std::make_unique<StateStorage>()),
                                                                          p_butcher(std::make_unique<ButcherTable>(stepper)),
                                                                          m_Nin(Nin),
                                                                          m_Nout(Nout),
                                                                          m_Npar(Npar){};

    Driver(int Nin, int Nout, int Npar) : p_states(std::make_unique<StateStorage>()),
                                          m_Nin(Nin),
                                          m_Nout(Nout),
                                          m_Npar(Npar){};

    ~Driver()
    {
    }

  public:
    int GetNin() const { return m_Nin; };
    int GetNout() const { return m_Nout; };
    int GetNpar() const { return m_Npar; };

    template <typename State>
    void GetState(State &u, int n) const { p_states->GetState(u, n); }

    int GetT() const { return p_states->GetT(); };

    template <typename Time>
    void PushBackTime(Time _time) { p_states->PushBackTime(_time); }

    template <typename State>
    void PushBackState(State state) { p_states->PushBackState(state); }

    double GetDt(int n) const { return p_states->GetDt(n); };
    double GetTime(int n) const { return p_states->GetTime(n); };

    // Public Member Functions
  public:
    template <typename State, typename Time>
    void Rhs(
        const State &u, State &dudt, const State &alphas, const Time &time)
    {
        p_aad_data->FeedWorkspace(u, alphas, time);

        p_aad_data->Forward();

        p_aad_data->ReadWorkspace(dudt);
    };
};

void delete_driver_handle(void *ptr)
{
    Driver *p_driver = static_cast<Driver *>(ptr);
    delete p_driver;
}

template <class Stepper>
void constructDriverButcherTableau(Driver &driver, Stepper stepper)
{
    // Construct Butcher Tableau that is needed for reverse pass
    driver.p_butcher = std::move(
        std::make_unique<ButcherTable>(stepper));
};

template <typename System>
void recordDriverRHSFunction(Driver &driver, System system)
{
    driver.p_aad_data = std::move(
        std::make_unique<AadData>(driver.GetNin(), driver.GetNpar(), system));
}

// Sets the derivatives of each component of the cost function w.r.t. to the ODE solution and w.r.t. the parameters
void setCostGradients(Driver &driver, std::vector<std::vector<double>> &lambda, std::vector<std::vector<double>> &mu)
{
    assert(lambda.size() == driver.GetNout());
    assert(lambda[0].size() == driver.GetNin());
    assert(mu.size() == driver.GetNout());
    assert(mu[0].size() == driver.GetNpar());

    driver.p_lambda = &lambda;
    driver.p_mu = &mu;

    // driver.p_mu = std::move(std::make_unique<std::vector<std::vector<double>>>(mu));
}

} // namespace lib

#endif
