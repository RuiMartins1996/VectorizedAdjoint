#ifndef _AADC_H_
#define _AADC_H_

#include <memory>
#include <chrono>
#include <map>
#include <aadc/idouble.h>
#include <aadc/ibool.h>
#include <aadc/iint.h>
#include <aadc/aadc_tools.h>
#include "aadc/aadc_ext.h"
#include <aadc/tools/runtime.h>
#include <aadc/tools/mmvector.h>
#include <aadc/aadc_thread_sync.h>
#include <aadc/aadc_stats.h>

typedef uint64_t AADCOpCounter;

extern "C" AADC_API AADCOpCounter CAAD_GetCurrentOpIndex();
extern "C" AADC_API void CAAD_ResetOptions();
extern "C" AADC_API void CAAD_SetOption(uint64_t option, uint64_t val);

extern "C" AADC_API void CAAD_WorkSpaceSetArray(int mmTypeSize, void* v, void* array_begin);

extern "C" AADC_API uint64_t CAAD_GetNumCheckPoints();
extern "C" AADC_API void CAAD_GetCheckPointCodeBlocks(uint64_t* check_point_blocks_indx);
extern "C" AADC_API uint64_t CAAD_GetNumCheckPointVars(uint64_t check_point_indx);
extern "C" AADC_API void CAAD_GetCheckPointVars(uint64_t check_point_indx, uint64_t* check_point_var_i);

extern "C" AADC_API int CAAD_iVarNumExtractPassive();

extern "C" AADC_API const char* CAAD_iVarGetExtractPassiveLocation(uint64_t loc_index);
extern "C" AADC_API uint64_t CAAD_iVarGetExtractPassiveLocationLine(uint64_t loc_index);

extern "C" AADC_API uint64_t CAAD_GetNumRegistredVars(int io_i);
extern "C" AADC_API void CAAD_GetRegisteredVar(int io_i, const char** var_id, bool* is_var_id_set, int64_t* i1, int64_t* i2, int64_t* i3, uint64_t* addr);

extern "C" AADC_API void CAAD_GetOpStats(uint64_t* fwd_op_stat_counter, uint64_t* rev_op_stat_counter);

extern "C" AADC_API const char* CAAD_Version();
extern "C" AADC_API const char* CAAD_Header_Checksums();

namespace aadc {

enum CompilerOption {
    AADC_MemOptimizations,
    AADC_RValueOptimizations,
    AADC_UseCompressedCode,
    AADC_NumCompressorThreads,
    AADC_CodeBufferSize,
    AADC_CodeCompressorFactor,
    AADC_CodeCompressorOrderFreeIndex,
    AADC_ReuseConstants,
    AADC_MaxNumVariables,
    AADC_InitInputDiff,
    AADC_ReleaseUncompressedCode,
    AADC_SingleAdjointPass,
    AADC_BreakOnActiveBoolConversion,
    AADC_InstrumentCode,
    AADC_DebugStopRecordingAt,
    AADC_LAST_OPTION
};

template<bool is_input>
class VarID {
public:
    VarID(
        const char* _var_id
        , int64_t _i1 = -1, int64_t _i2 = -1, int64_t _i3 = -1
    )
        : var_id(_var_id != 0 ? _var_id : "")
        , is_var_id_set(_var_id != 0)
        , i1(_i1), i2(_i2), i3(_i3)
    {}
    VarID(
        const std::string& _var_id
        , int64_t _i1 = -1, int64_t _i2 = -1, int64_t _i3 = -1
    )
        : var_id(_var_id)
        , is_var_id_set(true)
        , i1(_i1), i2(_i2), i3(_i3)
    {}
    VarID(
        int64_t _i1 = -1, int64_t _i2 = -1, int64_t _i3 = -1
    )
        : is_var_id_set(false)
        , i1(_i1), i2(_i2), i3(_i3)
    {}
    VarID(const VarID& other) = default;

    std::string var_id;
    bool is_var_id_set;
    int64_t i1, i2, i3;
};

template<bool is_input>
inline bool operator < (const VarID<is_input>& a, const VarID<is_input>& b) {
    bool eq = true;
    if (a.is_var_id_set < b.is_var_id_set) return true;
    eq = eq && a.is_var_id_set == b.is_var_id_set;
    if (eq && (a.i1 < b.i1)) return true;
    eq = eq && (a.i1 == b.i1);
    if (eq && (a.i2 < b.i2)) return true;
    eq = eq && (a.i2 == b.i2);
    if (eq && (a.i3 < b.i3)) return true;
    eq = eq && (a.i3 == b.i3);
    if (eq && a.is_var_id_set && b.is_var_id_set && (a.var_id < b.var_id)) return true;

    return false;
}

typedef VarID<true> ArgID;
typedef VarID<false> ResID;

template<typename mmType> class AADCFunctions;

class AADCFunctionArgRegistry {
    template<typename mmType> friend class AADCFunctions;
public:
    AADCArgument vArg(const std::string& var_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        ArgID look_for(var_id, i1, i2, i3);
        return vArg(look_for);
    }
    AADCArgument vArg(const ArgID& look_for) const {
        auto iter = input_var_reg.find(look_for);
        if (iter == input_var_reg.end()) throw std::runtime_error("AADC:Incorrect argument id requested");
        return iter->second;
    }
    AADCArgument dArg(const std::string& var_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        ArgID look_for(var_id, i1, i2, i3);
        return dArg(look_for);
    }
    AADCArgument dArg(const ArgID& look_for) const {
        auto iter = diff_var_reg.find(look_for);
        if (iter == diff_var_reg.end()) throw std::runtime_error("AADC:Incorrect diff argument id requested");
        return iter->second;
    }

    AADCResult vRes(const std::string& var_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        ResID look_for(var_id, i1, i2, i3);
        return vRes(look_for);
    }
    AADCResult vRes(const ResID& look_for) const {
        auto iter = output_var_reg.find(look_for);
        if (iter == output_var_reg.end()) throw std::runtime_error("AADC:Incorrect result id requested");
        return iter->second;
    }
private:
    std::map<ArgID, AADCArgument> input_var_reg;
    std::map<ArgID, AADCArgument> diff_var_reg;
    std::map<ResID, AADCResult> output_var_reg;
};


template<typename mmType>
class AADCWorkSpace {
    friend class AADCFunctions<mmType>;
public:
    AADCWorkSpace(
        const int64_t& work_array_size, const int64_t& stack_size
        , const std::shared_ptr<AADCFunctionArgRegistry>& arg_registry
    )
        : m_v(work_array_size)
        , m_d(work_array_size)
        , m_stack(stack_size+1)
        , m_arg_registry(arg_registry)
    {
      resetDiff();
    }
    AADCWorkSpace(
        const int64_t& work_array_size, const int64_t& stack_size
        , const std::shared_ptr<AADCFunctionArgRegistry>& arg_registry
        , const std::vector<uint64_t>& num_checkpoint_vars
    )
        : m_v(work_array_size)
        , m_d(work_array_size)
        , m_stack(stack_size+1)
        , m_check_points(num_checkpoint_vars.size())
        , m_arg_registry(arg_registry)
    {
      resetDiff();
      for (int i = 0; i < m_check_points.size(); ++i) {
          m_check_points[i].resize(num_checkpoint_vars[i]);
      }
    }
    AADCWorkSpace(const AADCWorkSpace& other)
        : m_v(other.m_v)
        , m_d(other.m_d)
        , m_stack(other.m_stack)
        , m_check_points(other.m_check_points)
        , m_rev_runs_after_fwd(other.m_rev_runs_after_fwd)
        , m_arg_registry(m_arg_registry)

    {}
    ~AADCWorkSpace() {}
public:
    void setArray(const mmType* array) {
        CAAD_WorkSpaceSetArray(sizeof(mmType), (void*)&(m_v[0]), (void*)array);
    }
    AADCWorkSpace& setVal(const AADCArgument& var_indx, const mmType& val) {
        m_v[var_indx.getIndex()] = val;
        return *this;
    }
    AADCWorkSpace& setVal(const AADCArgument& var_indx, const double val) {
        m_v[var_indx.getIndex()] = mmSetConst<mmType>(val);
        return *this;
    }
    AADCWorkSpace& setVal(const ArgID& var_id, const mmType& val) {
        return setVal(m_arg_registry->vArg(var_id), val);
    }
    AADCWorkSpace& setVal(const ArgID& var_id, const double& val) {
        return setVal(m_arg_registry->vArg(var_id), val);
    }
    AADCWorkSpace& setVal(const AADCScalarArgument& var_indx, const double val) {
        double* v((double*)(&m_v[0]));
        v[var_indx.getIndex()] = val;
        return *this;
    }

    template<typename IdxIter, typename ValIter>
    void setVal(IdxIter indx, ValIter val, ValIter val_end) {
        for (; val != val_end; ++indx, ++val) { setVal(*indx, *val); }
    }

    double& val(const AADCScalarArgument& var_indx) {
        double* v((double*)(&m_v[0]));
        return v[var_indx.getIndex()];
    }
    const double& val(const AADCScalarArgument& var_indx) const {
        double* v((double*)(&m_v[0]));
        return v[var_indx.getIndex()];
    }

    mmType& val(const AADCArgument& var_indx) {
        return m_v[var_indx.getIndex()];
    }
    const mmType& val(const AADCArgument& var_indx) const {
        return m_v[var_indx.getIndex()];
    }

    mmType& val(const ArgID& var_id) {
        return val(m_arg_registry->vArg(var_id));
    }
    const mmType& val(const ArgID& var_id) const {
        return val(m_arg_registry->vArg(var_id));
    }

    const double* valp(const AADCArgument& var_indx) const {
        return (double*)&(m_v[var_indx.getIndex()]);
    }
    double* valp(const AADCArgument& var_indx) {
        return (double*)&(m_v[var_indx.getIndex()]);
    }

    const double* valp(const ArgID& var_id) const {
        return valp(m_arg_registry->vArg(var_id));
    }
    double* valp(const ArgID& var_id) {
        return valp(m_arg_registry->vArg(var_id));
    }

    mmType& val(const AADCResult& var_indx) {
        return m_v[var_indx.getIndex()];
    }
    const mmType& val(const AADCResult& var_indx) const {
        return m_v[var_indx.getIndex()];
    }

    mmType& val(const ResID& var_id) {
        return val(m_arg_registry->vRes(var_id));
    }
    const mmType& val(const ResID& var_id) const {
        return val(m_arg_registry->vRes(var_id));
    }

    double* valp(const AADCResult& var_indx) {
        return (double*)&(m_v[var_indx.getIndex()]);
    }
    const double* valp(const AADCResult& var_indx) const {
        return (double*)&(m_v[var_indx.getIndex()]);
    }

    double* valp(const ResID& var_id) {
        return valp(m_arg_registry->vRes(var_id));
    }
    const double* valp(const ResID& var_id) const {
        return valp(m_arg_registry->vRes(var_id));
    }

    mmType& diff(const AADCArgument& var_indx) {
        return m_d[var_indx.getIndex()];
    }
    const mmType& diff(const AADCArgument& var_indx) const {
        return m_d[var_indx.getIndex()];
    }

    mmType& diff(const ArgID& var_id) {
        return diff(m_arg_registry->vArg(var_id));
    }
    const mmType& diff(const ArgID& var_id) const {
        return diff(m_arg_registry->vArg(var_id));
    }

    double* diffp(const AADCArgument& var_indx) {
        return (double*)&(m_d[var_indx.getIndex()]);
    }
    const double* diffp(const AADCArgument& var_indx) const {
        return (double*)&(m_d[var_indx.getIndex()]);
    }

    double* diffp(const ArgID& var_id) {
        return diffp(m_arg_registry->vArg(var_id));
    }
    const double* diffp(const ArgID& var_id) const {
        return diffp(m_arg_registry->vArg(var_id));
    }

    double& diff(const AADCScalarArgument& var_indx) {
        double* d((double*)(&m_d[0]));
        return d[var_indx.getIndex()];
    }
    const double& diff(const AADCScalarArgument& var_indx) const {
        double* d((double*)(&m_d[0]));
        return d[var_indx.getIndex()];
    }

	mmType& diff(const AADCResult& var_indx) {
        return m_d[var_indx.getIndex()];
    }
	const mmType& diff(const AADCResult& var_indx) const {
        return m_d[var_indx.getIndex()];
    }
	mmType& diff(const ResID& var_id) {
        return diff(m_arg_registry->vRes(var_id));
    }
	const mmType& diff(const ResID& var_id) const {
        return diff(m_arg_registry->vRes(var_id));
    }

    double* diffp(const AADCResult& var_indx) {
        return (double*)&(m_d[var_indx.getIndex()]);
    }
    const double* diffp(const AADCResult& var_indx) const {
        return (double*)&(m_d[var_indx.getIndex()]);
    }

    double* diffp(const ResID& var_id) {
        return diffp(m_arg_registry->vRes(var_id));
    }
    const double* diffp(const ResID& var_id) const {
        return diffp(m_arg_registry->vRes(var_id));
    }

    void setDiff(const AADCArgument& var_indx, const mmType& val) {
        m_d[var_indx.getIndex()] = val;
    }
    void setDiff(const AADCResult& var_indx, const mmType& val) {
        m_d[var_indx.getIndex()] = val;
    }
    void setDiff(const ResID& var_id, const mmType& val) {
        setDiff(m_arg_registry->vRes(var_id), val);
    }

    void setDiff(const AADCArgument& var_indx, const double val) {
        m_d[var_indx.getIndex()] = mmSetConst<mmType>(val);
    }
    void setDiff(const AADCResult& var_indx, const double val) {
        m_d[var_indx.getIndex()] = mmSetConst<mmType>(val);
    }
    void setDiff(const ResID& var_id, const double val) {
        setDiff(m_arg_registry->vRes(var_id), val);
    }

	// reset all adjoint variables
    void resetDiff() {
        std::fill(m_d.begin(), m_d.end(), mmZero<mmType>());
    }
private:
    mmVector<mmType> m_v;
    mmVector<mmType> m_d;
    mmVector<mmType> m_stack;
    std::vector<mmVector<mmType>> m_check_points;
    uint64_t m_rev_runs_after_fwd;
    std::shared_ptr<AADCFunctionArgRegistry> m_arg_registry;
};

template<typename mmType>
class AADCFunctions {
public:
    typedef void (* BinFunc)(mmType* v, mmType* d, const double* c, mmType* s);
    typedef std::vector<std::pair<CompilerOption, uint64_t>> AADCOptions;
public:
    AADCFunctions(const AADCOptions& options = AADCOptions()) 
        : m_options(options)
        , m_code_options(16)
    {
    }
    AADCFunctions(const AADCFunctions<mmType>& other) = delete;
    
    void setOption(CompilerOption option, uint64_t val) {
        m_options.push_back(std::pair<CompilerOption, uint64_t>(option, val));
    }
    void startRecording() {
        m_start_recording = std::chrono::high_resolution_clock::now();
        initOptions();
        if (CAAD_StartRecording() != 0) throw std::runtime_error("AADC license has expired");
        idouble::recording = true;
    }

    virtual void stopRecording() {
        m_arg_registry = std::make_shared<AADCFunctionArgRegistry>();

        for (int io_i = 0; io_i < 3; ++io_i) { // loop over inputs/outputs
            uint64_t num_reg_inputs(CAAD_GetNumRegistredVars(io_i));
            for (uint64_t i = 0; i < num_reg_inputs; ++i) {
                const char* var_id;
                bool is_var_id_set;
                int64_t i1, i2, i3;
                uint64_t addr;
                CAAD_GetRegisteredVar(io_i, &var_id, &is_var_id_set, &i1, &i2, &i3, &addr);

                if(io_i != 2) {
                    ArgID var_reg(var_id, i1, i2, i3);
                    m_arg_registry->input_var_reg[var_reg] = AADCArgument(addr);
                } else {
                    ResID var_reg(var_id, i1, i2, i3);
                    m_arg_registry->output_var_reg[var_reg] = AADCResult(addr);
                }
            }
        }
        
        CAAD_FinalizeRecording();

        m_work_array_size = CAAD_VDWorkArraySize();
        m_stack_size = CAAD_StackSize();
        m_const_data.resize(CAAD_ConstantsArraySize());
        double* c_vals = CAAD_ConstantsArray();
        std::copy(c_vals, c_vals + m_const_data.size(), m_const_data.begin());
        m_num_code_blocks = CAAD_GetNumCodeBlocks();
        m_func_forward.resize(m_num_code_blocks);
        m_func_reverse.resize(m_num_code_blocks);
        m_code_size_fwd = 0; m_code_size_rev = 0;
        for (uint64_t ci = 0; ci < m_num_code_blocks; ++ci) {
            uint64_t code_size_fwd = CAAD_GetCodeSize(true, ci);
            uint64_t code_size_rev = CAAD_GetCodeSize(false, ci);

            m_code_size_fwd += code_size_fwd; m_code_size_rev += code_size_rev;

            // TODO: go directly from ExtractCode to m_code
            std::vector<uint8_t > code(code_size_fwd);
            CAAD_ExtractCode(true, &(code[0]), ci);
            std::vector<uint8_t > rcode(code_size_rev);
            CAAD_ExtractCode(false, &(rcode[0]), ci);

//            std::cout << "CodeBlock : " << ci << " fwd code size : " << code_size_fwd 
//            << " rev code size : " << code_size_rev << std::endl;

            m_code.Add(&(m_func_forward[ci]), code);
            m_code.Add(&(m_func_reverse[ci]), rcode);
        }
        m_num_checkpoints = CAAD_GetNumCheckPoints();
        m_checkpoint_indx.resize(m_num_checkpoints);
        if (m_num_checkpoints > 0) CAAD_GetCheckPointCodeBlocks(&(m_checkpoint_indx.front()));
        m_op_code_func_start = CAAD_RecordingStartAt();
        m_op_code_func_end = CAAD_GetCurrentOpIndex()+1;

        m_checkpoint_vars.resize(m_num_checkpoints);m_num_checkpoint_vars.resize(m_num_checkpoints);
        m_total_num_checkpoint_vars = 0;
        for (int i = 0; i < m_num_checkpoints; ++i) {
            m_total_num_checkpoint_vars += m_num_checkpoint_vars[i] = CAAD_GetNumCheckPointVars(i);
            m_checkpoint_vars[i].resize(m_num_checkpoint_vars[i]);
            if (m_num_checkpoint_vars[i]> 0) CAAD_GetCheckPointVars(i, &(m_checkpoint_vars[i].front()));
        }
        m_ext_funcs_new = getRecordingConstStateExtFunctionRegistry();
        m_number_of_required_random_variables=CAAD_iVarNumberOfRequiredRandomVariables();
        m_stop_recording = std::chrono::high_resolution_clock::now();
        
        idouble::recording = false;
        CAAD_GetOpStats(&(fwd_op_stat_counter[0]), &(rev_op_stat_counter[0]));
    }

	virtual ~AADCFunctions() {}

    AADCFunctions& operator= (const AADCFunctions& other) = delete;
    
    std::shared_ptr<AADCWorkSpace<mmType> > createWorkSpace() const {
        if (m_num_checkpoints == 0) {
            return std::make_shared<AADCWorkSpace<mmType> >(
                m_work_array_size, m_stack_size, m_arg_registry
            );
        } else {
            return std::make_shared<AADCWorkSpace<mmType> >(
                m_work_array_size, m_stack_size, m_arg_registry, m_num_checkpoint_vars
            );
        }
    }
    void forward(
        AADCWorkSpace<mmType>& workspace, int first_cp = 0, int last_cp = -1,
        std::shared_ptr<AADCFunctionSynchronizer> thread_syncer = 0
    ) const {
        last_cp = last_cp == -1 ? int(m_num_checkpoints) : last_cp;
        if (m_num_checkpoints > 0) {
//            for (uint64_t cp = 0; cp <= m_num_checkpoints; ++cp) {
            for (uint64_t cp = first_cp; cp <= last_cp; ++cp) {
                uint64_t code_block_start = (cp == 0 ? 0 : m_checkpoint_indx[cp-1]);
                uint64_t code_block_end = (cp < m_num_checkpoints ? m_checkpoint_indx[cp] : getNumCodeBlocks());
                for (uint64_t ci = code_block_start; ci < code_block_end; ++ci) {
                    if (thread_syncer) thread_syncer->syncThreads();
                    ForwardFunc(ci)(
                        &(workspace.m_v.front()), &(workspace.m_d.front()),
                        &(m_const_data.front()), &(workspace.m_stack.front())
                    );
                }
                if (cp < m_num_checkpoints) {
                    saveCheckPointVars(cp, workspace);
                }
            }
        } else {
            for (uint64_t ci = 0; ci < getNumCodeBlocks(); ++ci) {
//                if (ext_func_call_iter < m_ext_func_code_block_index.size()) {
//                    if (m_ext_func_code_block_index[ext_func_call_iter] == ci) {
//                        m_ext_funcs[ext_func_call_iter]->forward(&(workspace.m_v.front()));
//                        ++ext_func_call_iter;
//                    }
//                }
                if (thread_syncer) thread_syncer->syncThreads();
                ForwardFunc(ci)(
                    &(workspace.m_v.front()), &(workspace.m_d.front()),
                    &(m_const_data.front()), &(workspace.m_stack.front())
                );
            }
        }
        workspace.m_rev_runs_after_fwd = 0;
    }
    void reverse(
        AADCWorkSpace<mmType>& workspace, int first_cp = 0, int last_cp = -1,
        std::shared_ptr<AADCFunctionSynchronizer> thread_syncer = 0
    ) {
        last_cp = last_cp == -1 ? int(m_num_checkpoints) : last_cp;
        if (m_num_checkpoints > 0) {
//            std::cout << "*******\n";
//            std::cout << "Num checkpoints " << m_num_checkpoints << std::endl;
//            for (int ci = 0; ci < m_checkpoint_indx.size(); ++ci) {
//                std::cout << " check at " <<  m_checkpoint_indx[ci];
//            }
//            std::cout << "\n*******\n";
            if (last_cp == m_num_checkpoints) {
                if (workspace.m_rev_runs_after_fwd > 0) {
    //                std::cout << "restore check point " << m_checkpoint_indx.size() - 1 << std::endl;
                    loadCheckPointVars(m_checkpoint_indx.size() - 1, workspace);
                    for (uint64_t ci = m_checkpoint_indx.back(); ci < getNumCodeBlocks(); ++ci) {
    //                    std::cout << "forward " << ci << std::endl;
                        if (thread_syncer) thread_syncer->syncThreads();
                        ForwardFunc(ci)(
                            &(workspace.m_v.front()), &(workspace.m_d.front()),
                            &(m_const_data.front()), &(workspace.m_stack.front())
                        );
                    }
                }
                for (int ci = int(getNumCodeBlocks())-1; ci >= int(m_checkpoint_indx.back()); --ci) {
    //                std::cout << "reverse " << ci << std::endl;
                    if (thread_syncer) thread_syncer->syncThreads();
                    ReverseFunc(ci)(
                        &(workspace.m_v.front()), &(workspace.m_d.front()),
                        &(m_const_data.front()), &(workspace.m_stack.front())
                    );
                }
                last_cp--;
            }
//            for (int checkp_i = int(m_checkpoint_indx.size())-1; checkp_i > 0; --checkp_i) {
            for (int checkp_i = last_cp; checkp_i >= first_cp; --checkp_i) {
                if (checkp_i > 0) {
    //                std::cout << "restore check point " << checkp_i << std::endl;
                    loadCheckPointVars(checkp_i, workspace);
    //                std::cout << "restore check point " << checkp_i-1 << std::endl;
                    loadCheckPointVars(checkp_i-1, workspace);
                    for (int ci = int(m_checkpoint_indx[checkp_i-1]); ci < int(m_checkpoint_indx[checkp_i]); ++ci) {
    //                    std::cout << "forward " << ci << std::endl;
                        if (thread_syncer) thread_syncer->syncThreads();
                        ForwardFunc(ci)(
                            &(workspace.m_v.front()), &(workspace.m_d.front()),
                            &(m_const_data.front()), &(workspace.m_stack.front())
                        );
                    }
                    for (int ci = int(m_checkpoint_indx[checkp_i])-1; ci >= int(m_checkpoint_indx[checkp_i-1]); --ci) {
    //                    std::cout << "reverse " << ci << std::endl;
                        if (thread_syncer) thread_syncer->syncThreads();
                        ReverseFunc(ci)(
                            &(workspace.m_v.front()), &(workspace.m_d.front()),
                            &(m_const_data.front()), &(workspace.m_stack.front())
                        );
                    }
                }
                else { // very first interval
                    for (int ci = int(m_checkpoint_indx[0])-1; ci >= 0; --ci) {
        //                std::cout << "reverse " << ci << std::endl;
                        loadCheckPointVars(0, workspace);
                        if (thread_syncer) thread_syncer->syncThreads();
                        ReverseFunc(ci)(
                            &(workspace.m_v.front()), &(workspace.m_d.front()),
                            &(m_const_data.front()), &(workspace.m_stack.front())
                        );
                    }
                }
            }
        } else {
//			int ext_func_call_iter(int(m_ext_func_code_block_index.size())-1);

            for (int64_t ci = getNumCodeBlocks()-1; ci >= 0; --ci) {
                if (thread_syncer) thread_syncer->syncThreads();
                ReverseFunc(ci)(
                    &(workspace.m_v.front()), &(workspace.m_d.front()),
                    &(m_const_data.front()), &(workspace.m_stack.front())
                );
//                if (ext_func_call_iter >= 0) {
//                    if (m_ext_func_code_block_index[ext_func_call_iter] == ci) {
//                        m_ext_funcs[ext_func_call_iter]->reverse(
//                            &(workspace.m_v.front()),
//                            &(workspace.m_d.front())
//                        );
//                        --ext_func_call_iter;
//                    }
//                }
            }
        }
        workspace.m_rev_runs_after_fwd++;
    }
    AADCOpCounter OpCodeStart() const {
        return m_op_code_func_start;
    }
    AADCOpCounter OpCodeEnd() const {
        return m_op_code_func_end;
    }
public:
    AADCArgument vArg(const std::string& var_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        return m_arg_registry->vArg(var_id, i1, i2, i3);
    }

    const std::map<ArgID, AADCArgument>& getInputArgumentIDs() const {
        return m_arg_registry->input_var_reg;
    }
    const std::map<ArgID, AADCArgument>& getDiffArgumentIDs() const {
        return m_arg_registry->diff_var_reg;
    }
    const std::map<ResID, AADCResult>& getResultIDs() const {
        return m_arg_registry->output_var_reg;
    }
public:
    uint64_t getCodeSizeFwd() const {
        return m_code_size_fwd;
    }
    uint64_t getCodeSizeRev() const {
        return m_code_size_rev;
    }
    uint64_t getWorkArraySize() const {
        return m_work_array_size;
    }
    uint64_t getStackSize() const {
        return m_stack_size;
    }
    uint64_t getConstDataSize() const {
        return m_const_data.size();
    }
    uint64_t getNumCodeBlocks() const {
        return m_num_code_blocks;
    }
    uint64_t getNumCheckPoints() const {
        return m_num_checkpoints;
    }
    uint64_t getNumCheckPointVars() const {
        return m_total_num_checkpoint_vars;
    }
    uint64_t getMemUse() const {
        return getCodeSizeFwd() + getCodeSizeRev() + getConstDataSize() * sizeof(double);
    }
    uint64_t getWorkSpaceMemUse() const {
        return (getWorkArraySize() * 2 + getStackSize() + getNumCheckPointVars()) * sizeof(mmType);
    }
    template<class ostream>
    void outStats(ostream& ostr, const char* func_name) {
        std::chrono::microseconds compile_time = std::chrono::duration_cast<std::chrono::microseconds>(m_stop_recording - m_start_recording);

        ostr << func_name << " Compile time      : " << compile_time.count() << std::endl;
        ostr << func_name << " Code size forward : " << getCodeSizeFwd() << std::endl;
        ostr << func_name << " Code size reverse : " << getCodeSizeRev() << std::endl;
        ostr << func_name << " Work array size   : " << getWorkArraySize() << std::endl;
        ostr << func_name << " Stack size        : " << getStackSize() << std::endl;
        ostr << func_name << " Const data size   : " << getConstDataSize() << std::endl;
        ostr << func_name << " CheckPoint size   : " << getNumCheckPointVars() << std::endl;
        ostr << func_name << " Func + WS mem use : " << getMemUse()+getWorkSpaceMemUse() << " bytes" << std::endl;
        ostr << func_name << " Num code blocks   : " << m_func_forward.size() << " " << m_func_reverse.size() << std::endl;

    }
    template<class ostream>
    void outOpStats(ostream& ostr, const char* func_name) {
        outFuncOpStats(ostr, func_name, "fwd", fwd_op_stat_counter);
        outFuncOpStats(ostr, func_name, "rev", rev_op_stat_counter);
    }
    template<class ostream>
    void outFuncOpStats(ostream& ostr, const char* func_name, const char* fwd_rev, uint64_t op_stat_counter[]) {
        const char** op_names(CAAD_GetOpsCounterNames());
        for (std::size_t opi = 0; opi < aadc::OpsCounter::AADC_LAST_OP; ++opi) {
            if (op_stat_counter[opi]) {
                ostr << func_name << " " << fwd_rev << " " << op_names[opi] << "\t: " << op_stat_counter[opi] << std::endl;
            }
        }
    }
    template<class ostream>
    void printPassiveExtractLocations(ostream& ostr, const char* func_name) {
        ostr << "Number active to passive conversions: " << CAAD_iVarNumExtractPassive() 
            << " while recording " << func_name
            << std::endl;
		for (int i = 0; i < CAAD_iVarNumExtractPassive(); ++i) {
			ostr << CAAD_iVarGetExtractPassiveLocation(i) << ":" << CAAD_iVarGetExtractPassiveLocationLine(i) << std::endl;
		}
    }
    uint64_t getNumberOfRequiredRandomVariables() {
        return m_number_of_required_random_variables;
    }
    const char* getAADCVersion() {
        return CAAD_Version();
    }
    const char* getAADCHeaderChecksums() {
        return CAAD_Header_Checksums();
    }
    template<class ostream>
    void versionFull(ostream& ostr) {
        ostr << getAADCVersion() << std::endl << getAADCHeaderChecksums() << std::endl;
    }

private:
    void initOptions() {
        CAAD_ResetOptions();
        for (auto i = m_options.begin(); i != m_options.end(); ++i) {
            CAAD_SetOption(i->first, i->second);
        }
        CAAD_Init(m_code_options, int(aadc::mmSize<mmType>()));
        getRecordingConstStateExtFunctionRegistry().resize(0);
    }
    virtual BinFunc ForwardFunc(uint64_t code_block) const {
        return m_func_forward[code_block];
    }
    virtual BinFunc ReverseFunc(uint64_t code_block) const {
        return m_func_reverse[code_block];
    }
    void saveCheckPointVars(uint64_t cp, AADCWorkSpace<mmType>& workspace) const {
        // Save checkpoint variables
        for (uint64_t cp_var_i = 0; cp_var_i < m_checkpoint_vars[cp].size(); ++cp_var_i) {
            uint64_t var_indx(m_checkpoint_vars[cp][cp_var_i]);
            workspace.m_check_points[cp][cp_var_i] = workspace.m_v[var_indx];
        }
    }
    void loadCheckPointVars(uint64_t cp, AADCWorkSpace<mmType>& workspace) const {
        for (uint64_t cp_var_i = 0; cp_var_i < m_checkpoint_vars[cp].size(); ++cp_var_i) {
            uint64_t var_indx(m_checkpoint_vars[cp][cp_var_i]);
            workspace.m_v[var_indx] = workspace.m_check_points[cp][cp_var_i];
        }
    }

protected:
    AADCOptions m_options;
    int m_code_options;
    std::vector<BinFunc> m_func_forward, m_func_reverse;
    std::shared_ptr<AADCFunctionArgRegistry> m_arg_registry;
    uint64_t m_code_size_fwd, m_code_size_rev;
    int64_t m_work_array_size;
    int64_t m_stack_size;
    uint64_t m_num_code_blocks;
    uint64_t m_total_num_checkpoint_vars;
    std::vector<uint64_t> m_checkpoint_indx;
    std::vector<double> m_const_data; // TODO: if recording multiple functions in one valuation flow, need to share const_data somehow
    JitRuntime m_code;

    uint64_t m_num_checkpoints;
    uint64_t m_number_of_required_random_variables;

    std::vector<std::vector<uint64_t>> m_checkpoint_vars;
    std::vector<uint64_t> m_num_checkpoint_vars;

//    std::vector<std::shared_ptr<ExtFunc> > m_ext_funcs;
//    std::vector<int> m_ext_func_code_block_index;

    std::vector<std::shared_ptr<ConstStateExtFunc> > m_ext_funcs_new;

    AADCOpCounter m_op_code_func_start, m_op_code_func_end;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start_recording, m_stop_recording;
    uint64_t fwd_op_stat_counter[aadc::AADC_LAST_OP];
    uint64_t rev_op_stat_counter[aadc::AADC_LAST_OP];
};
}; // napespace aadc


#endif
