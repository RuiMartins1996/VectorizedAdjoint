#ifndef AADDATA_NEW_HEADER // guard
#define AADDATA_NEW_HEADER

#include <aadc/aadc.h>
#include <aadc/aadc_debug.h>

// Boost library
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ublas = boost::numeric::ublas;
typedef __m256d mmType;

// TODO: This class needs to spawn several workspaces to allow parallelization

// template < class Stepper,class System, class State , class Time>
class AadData
{
  public:
    std::unique_ptr<aadc::AADCFunctions<mmType>> p_aad_funcs;
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> ws;

    std::vector<aadc::AADCArgument> input_args;
    std::vector<aadc::AADCArgument> input_parameters;
    std::unique_ptr<aadc::AADCArgument> input_time;
    std::vector<aadc::AADCResult> output_args;

    std::mutex m_mutex;

  public:
    //*CONSTRUCTOR
    template <class System>
    AadData(const int Nin, const int Nout, const int Npar, System system)
    {
        std::vector<aadc::AADCArgument> _input_args(Nin);
        input_args = _input_args;

        std::vector<aadc::AADCArgument> _input_parameters(Npar);
        input_parameters = _input_parameters;

        input_time = std::make_unique<aadc::AADCArgument>(0.0);

        std::vector<aadc::AADCResult> _output_args(Nout);
        output_args = _output_args;

        std::unique_ptr<aadc::AADCFunctions<mmType>> _p_aad_funcs =
            std::make_unique<aadc::AADCFunctions<mmType>>();
        p_aad_funcs = std::move(_p_aad_funcs);

        ws = Record(p_aad_funcs, system);
    };

    // Destructor
    ~AadData() = default;

    // Delete copy constructors
    AadData(AadData const &) = delete;
    AadData &operator=(AadData const &) = delete;

    // Move assignment
    AadData &operator=(AadData &&other)
    {
        if (this != &other) {
            p_aad_funcs = std::move(other.p_aad_funcs);
            ws = std::move(other.ws);
            input_args = std::move(other.input_args);
            input_parameters = std::move(other.input_parameters);
            input_time = std::move(other.input_time);
            output_args = std::move(other.output_args);
        }
        return *this;
    };
    // Move Constructor
    AadData(AadData &&other) noexcept
    {
        p_aad_funcs = std::move(other.p_aad_funcs);
        ws = std::move(other.ws);
        input_args = std::move(other.input_args);
        input_parameters = std::move(other.input_parameters);
        input_time = std::move(other.input_time);
        output_args = std::move(other.output_args);
    };

    //* ACCESSORS
    static constexpr int n_d_in_mmtype()
    {
        return static_cast<int>(sizeof(mmType) / sizeof(double));
    };

    int GetInputArgsSize() const { return input_args.size(); };
    int GetInputParametersSize() const { return input_parameters.size(); };
    int GetOutputArgsSize() const { return output_args.size(); };

    std::vector<aadc::AADCArgument> GetInputArgs() const { return input_args; }
    std::vector<aadc::AADCArgument> GetInputPar() const { return input_parameters; }
    aadc::AADCArgument GetInputTime() const { return *input_time; }
    std::vector<aadc::AADCResult> GetOutputArgs() const { return output_args; }

    //*FORWARD:
    void Forward() { p_aad_funcs->forward(*(ws)); }
    //* REVERSE
    void Reverse() { p_aad_funcs->reverse(*(ws)); }

    //*RECORD : Calls the record function from the AADC library on the
    // user-defined rhs function
    template <class System>
    std::shared_ptr<aadc::AADCWorkSpace<mmType>> Record(
        std::unique_ptr<aadc::AADCFunctions<mmType>> &p_aad_funcs,
        System system)
    {
        int Nin = input_args.size();
        int Npar = input_parameters.size();
        int Nout = output_args.size();

        //! No need for initialization
        // Declare state vector as idouble
        std::vector<idouble> iy(Nin);
        for (int i = 0; i < Nin; i++)
            iy[i] = idouble(i);
        // Declare parameter vector as idouble
        std::vector<idouble> ialphas(Npar, idouble(1.0));
        // Declare output state vector as idouble
        std::vector<idouble> iF(iy.size());
        // Declare time as idouble
        idouble it;

        //* Record function
        p_aad_funcs->startRecording();

        // Mark state vector as inputs
        for (int i = 0; i < Nin; i++) {
            input_args[i] = iy[i].markAsInput();
        } // vector input

        // Mark parameter vector as inputs
        for (int i = 0; i < Npar; i++) {
            input_parameters[i] = ialphas[i].markAsInput();
        } // vector input

        // Mark time as input
        *input_time = it.markAsInput();

        // Run specified function
        system(iy, iF, ialphas, it);

        // Mark function result as output
        for (int i = 0; i < Nout; i++) {
            output_args[i] = iF[i].markAsOutput();
        }
        p_aad_funcs->stopRecording();

        return p_aad_funcs->createWorkSpace();
    };

    //*FEEDWORKSPACE: Takes a state vector, parameters and time and inserts it
    // on the AADC workspace
    template <typename State, typename Value>
    void FeedWorkspace(
        // Inputs
        const State &y, const State &alpha, const Value &t)
    {

        // Copy the same values to every AVX entry
        int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));

        // Feed state vector
        for (int i = 0; i < y.size(); i++) {
            for (int j = 0; j < sizeOfBatch; j++) {
                ws->val(input_args[i])[j] = y[i];
            }
        }

        // Feed parameters vector
        for (int i = 0; i < alpha.size(); i++) {
            for (int j = 0; j < sizeOfBatch; j++) {
                ws->val(input_parameters[i])[j] = alpha[i];
            }
        };

        // Feed time
        for (int j = 0; j < sizeOfBatch; j++) {
            ws->val(*input_time)[j] = t;
        }
    }

    template <class State>
    void ReadWorkspace(State &dudt) const
    {
        // Read the workspace output
        for (int i = 0; i < output_args.size(); i++)
            dudt[i] = ws->val(output_args[i])[0];
    };

    ublas::matrix<double> ReadWorkspace() const
    {
        int N = output_args.size();
        ublas::matrix<double> fy(N, 1); //!(should be(N,Ndiv))

        // Read the workspace output
        for (int i = 0; i < N; i++)
            //! for (int j = 0; j < Ndiv; j++) (REMOVED THIS)
            fy(i, 0) = ws->val(output_args[i])[0];
        return fy;
    }

    template <typename State, typename Value, typename Matrix>
    void FuncDerivatives( // Inputs
        const State &y, const State &alpha, const Value &time,
        // Outputs
        Matrix &dFdw, Matrix &dFdalpha)
    {
        std::lock_guard<std::mutex> guard(m_mutex);
        int Nin = this->input_args.size();
        int Npar = this->input_parameters.size();
        int Nout = this->output_args.size();

        // Feed workspace and run function forward
        FeedWorkspace(y, alpha, time);

        p_aad_funcs->forward(*(this->ws));

        // Reset the diff storage (HAVE TO INITIALIZE THEM TO ZERO BEFORE
        // CHANGING!)
        this->ws->resetDiff();

        // Make sure all are set to zero (resetDiff should do this!)
        for (int i = 0; i < Nout; i++) {
            this->ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(0.0);
        }

        // Obtain derivatives with respect to parameter vector entries and store
        for (int i = 0; i < Nout; i++) {
            double value;
            // Seed to compute adjoints with respect to output i
            this->ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(1.0);

            // Call AADC to get adjoints
            p_aad_funcs->reverse(*(this->ws));

            for (int j = 0; j < Npar; j++) {
                value = this->ws->diff(input_parameters[j])[0];
                if (isnan(value)) {
                    value = 0.0;
                };
                dFdalpha(i, j) = value;
            }

            // Reset seed
            this->ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(0.0);
        }

        // Obtain derivatives with respect to state vector entries and store
        for (int i = 0; i < Nout; i++) {
            double value = 0.0;
            // Seed to compute adjoints with respect to output i
            this->ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(1.0);
            // Call AADC to get adjoints
            p_aad_funcs->reverse(*(this->ws));

            for (int j = 0; j < Nin; j++) {
                value = this->ws->diff(input_args[j])[0];
                if (isnan(value)) {
                    value = 0.0;
                };
                dFdw(i, j) = value;
            }

            // Reset seed
            this->ws->diff(output_args[i]) = aadc::mmSetConst<mmType>(0.0);
        }
    }

    template <class State, class Value>
    void VectorTimesJacobians(
        // Inputs
        const aadc::mmVector<mmType> &w, const State &y, const State &alpha,
        const Value &time,
        // Outputs
        aadc::mmVector<mmType> &wtimesJx,
        aadc::mmVector<mmType> &wtimesJalpha)
    {

        std::lock_guard<std::mutex> guard(m_mutex);

        // std::cout<<"Worker "<< std::this_thread::get_id() <<"computing Jacobians"<<std::endl;

        // Feed workspace and run function forward
        FeedWorkspace(y, alpha, time);
        p_aad_funcs->forward(*(this->ws));

        this->ws->resetDiff();

        // Initialize seed to values of w
        for (int i = 0; i < output_args.size(); i++) {
            this->ws->diff(this->output_args[i]) = w[i];
        }

        // Call AADC to get adjoints
        p_aad_funcs->reverse(*(this->ws));

        // Adjoints for the input arguments
        for (int i = 0; i < input_args.size(); i++) {
            wtimesJx[i] = this->ws->diff(this->input_args[i]);
        }

        // Adjoints for the input parameters
        for (int k = 0; k < input_parameters.size(); k++) {
            wtimesJalpha[k] = this->ws->diff(this->input_parameters[k]);
        }

        // std::cout<<"Worker "<< std::this_thread::get_id() <<"finished"<< std::endl;
    };

    template <class State, class Value>
    void VectorTimesJacobians(
        // Inputs
        const State &w, const State &y, const State &alpha, const Value &time,
        // Outputs
        State &wtimesJx, State &wtimesJalpha)
    {

        std::lock_guard<std::mutex> guard(m_mutex);
        // Feed workspace and run function forward
        FeedWorkspace(y, alpha, time);
        p_aad_funcs->forward(*(this->ws));

        this->ws->resetDiff();

        // Initialize seed to values of w
        for (int i = 0; i < output_args.size(); i++) {
            this->ws->diff(this->output_args[i]) =
                aadc::mmSetConst<mmType>(w[i]);
        }

        // Call AADC to get adjoints
        p_aad_funcs->reverse(*(this->ws));

        // Adjoints for the input arguments
        for (int j = 0; j < input_args.size(); j++) {
            wtimesJx[j] = this->ws->diff(this->input_args[j])[0];
            if (isnan(wtimesJx[j])) {
                wtimesJx[j] = 0.0;
            }
        }

        // Adjoints for the input parameters
        for (int k = 0; k < input_parameters.size(); k++) {
            wtimesJalpha[k] = this->ws->diff(this->input_parameters[k])[0];
            if (isnan(wtimesJalpha[k])) {
                wtimesJalpha[k] = 0.0;
            }
        }

        //! Correct for Nans?!
    };
};

template <class System>
std::unique_ptr<AadData> AadDataFactory(const int Nin, const int Nout, const int Npar,
                                        System system)
{
    // Implicit move operation into the variable that stores the result.
    // AadData aad(Nin, Npar, stepper, system);
    return std::make_unique<AadData>(Nin, Nout, Npar, system);
};

#endif
