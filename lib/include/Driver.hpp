#ifndef ODEDRIVER_NEW_HPP // guard
#define ODEDRIVER_NEW_HPP

#include "AadData.hpp"
#include "ButcherTable.hpp"
#include "StateStorage.hpp"
#include "utilities.hpp"

// Boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

class Driver
{
  public:
    std::unique_ptr<AadData> p_aad_data;
    std::unique_ptr<StateStorage> p_states;
    std::unique_ptr<ButcherTable> p_butcher;

    int m_Nin;
    int m_Npar;

  public:
    template <class Stepper, class System, class State>
    Driver(Stepper stepper, System system, int Nin, int Nout, int Npar, State &u0) : p_aad_data(std::make_unique<AadData>(Nin, Nout, Npar, system)),
                                                                                     p_states(std::make_unique<StateStorage>()),
                                                                                     p_butcher(std::make_unique<ButcherTable>(stepper, u0)),
                                                                                     m_Nin(Nin),
                                                                                     m_Npar(Npar){};

    ~Driver() {}

  public:
    int GetNin() const { return m_Nin; };
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

#endif
