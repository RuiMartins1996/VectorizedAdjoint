


template<class mmType>
class Jacobian {
public:
    Jacobian(size_t num_diff, size_t num_out)
        : m_num_diff(num_diff)
        , m_num_out(num_out)
    {}

private:
    size_t m_num_diff, m_num_out;
    mmVector<mmType> m_data;
};

template<class mmType>
class InputAdjoints {
public:
    InputAdjoints(size_t num_diff)
        : m_num_diff(num_diff)
    {}

    const mmType& operator ()(const aadc::AADCArgument& arg) {
        m_data[arg.getDiffIndex()];
    }

private:
    size_t m_num_diff;
    mmVector<mmType> m_data;
};


template<class mmType, bool accumulate>
void calcJacobian(
    vector<InputAdjoints<mmType> >& jac, WorkSpace<mmType>& ws
    , const AADCFunction<mmType>& func
    , const std::vector<aadc::AADCResult>& outputs
) {
    jac.resize(output.size());

    func.forward(ws);

    for (int i = 0; i < output.size(); ++i) {
        ws.diff(outputs[i]) = mmZero<mmType>();
        jac[i].resize(func.getNumDiffVars());
    }

    for (int i = 0; i < output.size(); ++i) {
        ws.diff(outputs[i]) = mmSetConst<mmType>(1.0);
        func.reverse(ws);
        for (int j = 0; j < func.getNumDiffVars(); ++j) {
            if (accumulate) {
                jac[i][j] = mmAdd(jac[i][j], ws.diff(func.getDiffVar(j)));
            } else {
                jac[i][j] = ws.diff(func.getDiffVar(j));
            }
        }
    }
}