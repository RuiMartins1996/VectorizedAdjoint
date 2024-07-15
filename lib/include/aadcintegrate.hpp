#ifndef AADC_INTEGRATE_HPP_INCLUDED
#define AADC_INTEGRATE_HPP_INCLUDED

#include <aadc/aadc.h>
#include <aadc/aadc_debug.h>

#include "runge_kutta.hpp"
#include "SystemFunctor.hpp"

namespace odeint = boost::numeric::odeint;
typedef __m256d mmType;

class AadDataForRK {
   public:
    aadc::AADCFunctions<mmType> aad_funcs;
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> ws;

    std::vector<aadc::AADCArgument> input_args;
    std::vector<aadc::AADCArgument> input_parameters;
    std::vector<aadc::AADCResult> output_args;

   public:
    //*CONSTRUCTOR
    template <class System>
    AadDataForRK(const int Nin, const int Npar, System system, double ti, double tf, double dt) {
        std::vector<aadc::AADCArgument> _input_args(Nin);
        input_args = _input_args;

        std::vector<aadc::AADCArgument> _input_parameters(Npar);
        input_parameters = _input_parameters;

        std::vector<aadc::AADCResult> _output_args(Nin);
        output_args = _output_args;

        ws = Record(aad_funcs, system, ti, tf, dt);
    };

    template <class System>
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> Record(
        aadc::AADCFunctions<mmType> &aad_funcs, System system, double ti, double tf, double dt) {
        typedef boost::numeric::odeint::runge_kutta_cash_karp54<
            std::vector<idouble>>
            stepper_type;

        stepper_type stepper;

        int Nin = input_args.size();
        int Npar = input_parameters.size();
        int Nout = output_args.size();

        // Declare state vector as idouble
        std::vector<idouble> iy(Nin);
        // Declare parameter vector as idouble
        std::vector<idouble> ialphas(Npar);
        // Declare output state vector as idouble
        std::vector<idouble> iF(iy.size());

        //* Record function
        aad_funcs.startRecording();

        // Mark state vector as inputs
        for (int i = 0; i < Nin; i++) {
            input_args[i] = iy[i].markAsInput();
        }  // vector input

        // Mark parameter vector as inputs
        for (int i = 0; i < Npar; i++) {
            input_parameters[i] = ialphas[i].markAsInput();
        }  // vector input

        auto system_runge_kutta_step =
            SystemFunctor<System, std::vector<idouble>, idouble>(system,
                                                                 ialphas);

        stepper_type st;

        integrate_const(st, system_runge_kutta_step, iy, ti, tf, dt);

        // Mark function result as output
        for (int i = 0; i < Nout; i++) {
            iF[i] = iy[i];
            output_args[i] = iF[i].markAsOutput();
        }
        aad_funcs.stopRecording();

        return aad_funcs.createWorkSpace();
    };
    //* ACCESSORS
    int GetInputArgsSize() const { return input_args.size(); };
    int GetInputParametersSize() const { return input_parameters.size(); };
    int GetOutputArgsSize() const { return output_args.size(); };

    std::vector<aadc::AADCArgument> GetInputArgs() { return input_args; }
    std::vector<aadc::AADCArgument> GetInputPar() { return input_parameters; }
    std::vector<aadc::AADCResult> GetOutputArgs() { return output_args; }

    //*FORWARD:
    void Forward() { aad_funcs.forward(*(ws)); }
    //* REVERSE
    void Reverse() { aad_funcs.reverse(*(ws)); }

    //*FEEDWORKSPACE: Takes a state vector, parameters and time and inserts it
    // on the AADC workspace
    template <typename State>
    void FeedWorkspace(
        // Inputs
        const State &y, const State &alpha) {
        // Copy the same values to every AVX entry

        //! for (int j = 0; j < Ndiv; j++)(REMOVED THIS!)

        // Feed state vector
        for (int i = 0; i < y.size(); i++) {
            ws->val(input_args[i])[0] = y[i];
        }

        // Feed parameters vector
        for (int i = 0; i < alpha.size(); i++) {
            ws->val(input_parameters[i])[0] = alpha[i];
        };
    }

    template <typename State>
    void Run(const State &y, State &dudt, const State &alpha) {
        FeedWorkspace(y, alpha);
        Forward();
        ReadWorkspace(dudt);
    }

    template <class State>
    void ReadWorkspace(State &dudt) const {
        // Read the workspace output
        for (int i = 0; i < output_args.size(); i++)
            dudt[i] = ws->val(output_args[i])[0];
    };
    template <class State>
    void VectorTimesJacobians(
        // Inputs
        const State &w, const State &y, const State &alpha,
        // Outputs
        State &wtimesJx, State &wtimesJalpha) {
        int Nin = input_args.size();
        int Npar = input_parameters.size();
        int Nout = output_args.size();

        // Feed workspace and run function forward
        FeedWorkspace(y, alpha);
        aad_funcs.forward(*(this->ws));

        this->ws->resetDiff();

        // Initialize seed to values of w
        for (int i = 0; i < Nout; i++) {
            this->ws->diff(this->output_args[i]) =
                aadc::mmSetConst<mmType>(w[i]);
        }

        // Call AADC to get adjoints
        aad_funcs.reverse(*(this->ws));

        // Adjoints for the input arguments
        for (int j = 0; j < Nin; j++) {
            wtimesJx[j] = this->ws->diff(this->input_args[j])[0];
        }

        // Adjoints for the input parameters
        for (int k = 0; k < Npar; k++) {
            wtimesJalpha[k] = this->ws->diff(this->input_parameters[k])[0];
        }
    };
};

#endif