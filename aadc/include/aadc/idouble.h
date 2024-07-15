#ifndef CAAD_IDOUBLE_H
#define CAAD_IDOUBLE_H

#include <stdint.h>
#include <cmath>
#include <math.h>
#include <limits>
#include <ostream>

// #define AADC_QL
#ifndef AADC_PASSIVE
#define AADC_IDOUBLE_ACTIVE
#endif

#include <aadc/api.h>

//
// define AADC_CHECK_VALUES to make debugger stop if there are any nan or inf values during
// calculations/recording
//
// #define AADC_CHECK_VALUES

#ifdef AADC_CHECK_VALUES
    #ifdef _MSC_VER
    #define AADC_CHECK_VALUE(x) if (!isfinite(x)) __debugbreak();
    #else
    #include <signal.h>
    #define AADC_CHECK_VALUE(x) if (!isfinite(x)) raise(SIGTRAP);
    #endif
#else
    #define AADC_CHECK_VALUE(x)
#endif

#define FM_DOUBLE_LOG2OFE  1.4426950408889634074
#define FM_DOUBLE_LOGEOF2  6.9314718055994530942e-1
#define FM_DOUBLE_LOG2OF10 3.32192809488736234789

extern "C" AADC_API void CAAD_Init(int cf, int numv);
extern "C" AADC_API void CAAD_iVarConstructor(double*);
extern "C" AADC_API void CAAD_iVarConstructorUninitialized(double*);
extern "C" AADC_API void CAAD_iVarConstructorConstant(double*, const double*);

extern "C" AADC_API void CAAD_iVarDestructor(const double*);

extern "C" AADC_API void CAAD_iInitVarWithSpecificAddress(const double* v, bool scalar, uint64_t given_addr);
extern "C" AADC_API uint64_t CAAD_iReserveContinuousAddressSpace(bool scalar, uint64_t size);

extern "C" AADC_API void CAAD_iVarMarkAsInput(const double*, const char* str_id, int64_t i1, int64_t i2, int64_t i3);
extern "C" AADC_API void CAAD_iVarMarkAsDiff(const double*, const char* str_id, int64_t i1, int64_t i2, int64_t i3);
extern "C" AADC_API void CAAD_iVarMarkAsInputNoDiff(const double*, const char* str_id, int64_t i1, int64_t i2, int64_t i3);
extern "C" AADC_API void CAAD_iVarMarkAsOutput(const double*, const char* str_id, int64_t i1, int64_t i2, int64_t i3);

extern "C" AADC_API void CAAD_iVarMarkAsScalarInput(const double* v);

extern "C" AADC_API void CAAD_iVarForceVariable(const double* v);
extern "C" AADC_API void CAAD_iVarForceDiff(const double* v);

extern "C" AADC_API void CAAD_iVarMarkAsArrayInput(const double* v, uint32_t indx);
extern "C" AADC_API void CAAD_iVarMarkAsArrayInputNoDiff(const double* v, uint32_t indx);

extern "C" AADC_API void CAAD_iVarSaveToStack(const double*);
extern "C" AADC_API double CAAD_iVarExtractPassive(const double*, const char*, uint64_t line);


extern "C" AADC_API void CAAD_iVarAssign(const double* in, double* out, bool rvalue);

extern "C" AADC_API void CAAD_iVarAddAssign(const double* in, double* out);
extern "C" AADC_API void CAAD_iVarAddAssignR(const double* in, double* out);
extern "C" AADC_API void CAAD_iVarSubAssign(const double* in, double* out);
extern "C" AADC_API void CAAD_iVarAdd(const double* in1, const double* in2, double* out, const bool rv1, const bool rv2);
//extern "C" AADC_API void CAAD_iVarAdd(const double* in1, const double* in2, double* out);
extern "C" AADC_API void CAAD_iVarSub(const double* in1, const double* in2, double* out, const bool rv1, const bool rv2);
//extern "C" AADC_API void CAAD_iVarSub(const double* in1, const double* in2, double* out);
extern "C" AADC_API void CAAD_iVarMult(const double* in1, const double* in2, double* out, const bool rv1, const bool rv2);
//extern "C" AADC_API void CAAD_iVarMult(const double* in1, const double* in2, double* out);
extern "C" AADC_API void CAAD_iVarDiv(const double* in1, const double* in2, double* out);
extern "C" AADC_API void CAAD_iVarCDFNorm(const double* in, double* out, double* out_pdf);

extern "C" AADC_API void CAAD_iVarFmod(const double* in1, const double* in2, double* out);

extern "C" AADC_API void CAAD_iVarFuncCall1Arg(const char* name, uint64_t func_adr, const double* in, double* out);

extern "C" AADC_API void CAAD_CodeFlush();
extern "C" AADC_API int64_t CAAD_StartRecording();
extern "C" AADC_API uint64_t CAAD_RecordingStartAt();

extern "C" AADC_API void CAAD_CheckPoint();
extern "C" AADC_API void CAAD_DoNotZeroInitAllCurrentDiffVariables();

extern "C" AADC_API bool CAAD_Output_Cpp(const char* path_output, const char* name, uint64_t code_block);
extern "C" AADC_API bool CAAD_Output_Asm(const char* path_output, const char* name, uint64_t code_block);

extern "C" AADC_API void CAAD_FinalizeRecording();

extern "C" AADC_API uint64_t CAAD_GetNumCodeBlocks();

extern "C" AADC_API uint64_t CAAD_GetAVXSize();

extern "C" AADC_API void CAAD_ExtractCode(bool forward, uint8_t* codep, uint64_t code_block);
extern "C" AADC_API uint64_t CAAD_GetCodeSize(bool forward, uint64_t code_block);

extern "C" AADC_API uint64_t CAAD_ConstantsArraySize();
extern "C" AADC_API double* CAAD_ConstantsArray();

extern "C" AADC_API uint64_t CAAD_VDWorkArraySize();

extern "C" AADC_API uint64_t CAAD_StackSize();

extern "C" AADC_API uint64_t CAAD_iVarIndex(const double*);

extern "C" AADC_API void CAAD_Assert(const char* msg, const double* val);

extern "C" AADC_API void CAAD_LoopPulse(const uint64_t loop_pulse);
extern "C" AADC_API void CAAD_StopPulses();

extern "C" AADC_API void CAAD_Comment(const char* comment);

extern "C" AADC_API bool CAAD_iVarIsRandom(const double *);
extern "C" AADC_API bool CAAD_iVarIsDiff(const double *);

extern "C" AADC_API void CAAD_iVarMarkVariableAsRandomInput(const double* v);
extern "C" AADC_API uint64_t CAAD_iVarNumberOfRequiredRandomVariables();

extern "C" AADC_API void AADCDebugInstrumentationOn();
extern "C" AADC_API void AADCDebugInstrumentationOff();

extern "C" AADC_API void AADCDebugPrintOn();
extern "C" AADC_API void AADCDebugPrintOff();

class idouble;

inline idouble operator - (const idouble& a, const idouble& b) noexcept;

#ifdef _WIN64
#define AADC_LIKELY(x) (x)
#define AADC_UNLIKELY(x) (x)
#else
#define AADC_LIKELY(x) __builtin_expect(!!(x), 1)
#define AADC_UNLIKELY(x) __builtin_expect(!!(x), 0)
#endif
namespace aadc {

class AADCArgument {
public:
    explicit AADCArgument(uint64_t index)
        : m_index(index)
    {}
    AADCArgument() = default;
    AADCArgument(const AADCArgument& other) = default;
    AADCArgument& operator = (const AADCArgument& other) = default;
    uint64_t getIndex() const {
        return m_index;
    }    
private:
    uint64_t m_index;
};

class AADCScalarArgument {
public:
    explicit AADCScalarArgument(uint64_t index)
        : m_index(index)
    {}
    AADCScalarArgument() = default;
    AADCScalarArgument(const AADCScalarArgument& other) = default;
    AADCScalarArgument& operator = (const AADCScalarArgument& other) = default;
    uint64_t getIndex() const {
        return m_index;
    }    
private:
    uint64_t m_index;
};

class AADCVectorArgument {
public:
    AADCVectorArgument(uint64_t index, uint64_t len)
        : m_index(index)
        , m_len(len)
    {}
    AADCVectorArgument() = default;
    AADCVectorArgument(const AADCVectorArgument& other) = default;
    AADCVectorArgument& operator = (const AADCVectorArgument& other) = default;
    uint64_t getIndex() const {
        return m_index;
    }
    uint64_t size() const {
        return m_len;
    }
private:
    uint64_t m_index, m_len;
};

class AADCScalarVectorArgument {
public:
    AADCScalarVectorArgument(uint64_t index, uint64_t len)
        : m_index(index)
        , m_len(len)
    {}
    AADCScalarVectorArgument() = default;
    AADCScalarVectorArgument(const AADCScalarVectorArgument& other) = default;
    AADCScalarVectorArgument& operator = (const AADCScalarVectorArgument& other) = default;
    uint64_t getIndex() const {
        return m_index;
    }
    uint64_t size() const {
        return m_len;
    }
private:
    uint64_t m_index, m_len;
};

class AADCResult {
public:
    AADCResult(uint64_t index)
        : m_index(index)
    {}
    AADCResult() = default;
    AADCResult(const AADCResult& other) = default;
    AADCResult& operator = (const AADCResult& other) = default;

    uint64_t getIndex() const {
        return m_index;
    }

private:
    uint64_t m_index;
};

}; // namespace aadc

inline std::ostream& operator<<(std::ostream& os, const aadc::AADCArgument& a)
{
    os << "[Input, " << a.getIndex() << "]";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const aadc::AADCResult& a)
{
    os << "[Result, " << a.getIndex() << "]";
    return os;
}

class idouble {
friend inline std::ostream& operator<<(std::ostream&, const idouble&);
public:
AADC_API static bool recording;
AADC_API static int64_t DebugOpCounter;
AADC_API static int64_t DebugOpCounterStart;
AADC_API static int64_t DebugOpCounterEnd;
public:
    idouble(const double& c) noexcept : val(c) {
        AADC_CHECK_VALUE(val)
//        if(AADC_UNLIKELY(recording)) CAAD_iVarConstructorConstant(&val, &c);
    }
    idouble(idouble&& other)  noexcept
        : val(other.val)
    {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(other.val), &val, true);
#endif
        AADC_CHECK_VALUE(val);
    }
    idouble()  noexcept : val(0.0) { 
//        if(AADC_UNLIKELY(recording)) CAAD_iVarConstructorUninitialized(&val); 
    }
    idouble(const idouble& other) 
        : val(other.val) 
    { 
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&other.val, &val, false);
#endif
        AADC_CHECK_VALUE(val);
    }
    ~idouble()  noexcept {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarDestructor(&val);
#endif
    }

    idouble& operator = (const idouble& other) noexcept {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(other.val), &val, false);
#endif
        val = other.val;
        AADC_CHECK_VALUE(val);
        return *this;
    }
    idouble& operator = (const idouble&& other) noexcept {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(other.val), &val, true);
#endif
        val = other.val;
        AADC_CHECK_VALUE(val);
        return *this;
    }
    idouble& operator += (const idouble& other) noexcept {
        idouble tmp;
        // TODO: Implement += more efficiently
        tmp.val = val + other.val;
#ifdef AADC_IDOUBLE_ACTIVE
        if(AADC_UNLIKELY(recording)) CAAD_iVarAdd(&val, &other.val, &tmp.val, false, false);
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(tmp.val), &val, false);
//        if(AADC_UNLIKELY(recording)) CAAD_iVarAddAssign(&(other.val), &val);
#endif
        val = tmp.val;
        AADC_CHECK_VALUE(val);
        return *this;
    }
    idouble& operator += (const idouble&& other) noexcept {
        idouble tmp;
        // TODO: Implement += more efficiently
        tmp.val = val + other.val;
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarAdd(&val, &other.val, &tmp.val, false, false);
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(tmp.val), &val, false);
//        if(AADC_UNLIKELY(recording)) CAAD_iVarAddAssignR(&(other.val), &val);
#endif
        val = tmp.val;
        AADC_CHECK_VALUE(val);
        return *this;
    }
    idouble& operator -= (const idouble& other) noexcept {
        idouble tmp;
        // TODO: Implement -= more efficiently
        tmp.val = val - other.val;
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarSub(&val, &other.val, &tmp.val, false, false);
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(tmp.val), &val, false);
//        if(AADC_UNLIKELY(recording)) CAAD_iVarSubAssign(&(other.val), &val);
#endif
        val = tmp.val;
        AADC_CHECK_VALUE(val);
        return *this;
    }
    idouble& operator *= (const idouble& other) noexcept {
        idouble tmp;
        // TODO: Implement *= more efficiently
        tmp.val = val * other.val;
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarMult(&val, &other.val, &tmp.val, false, false);
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(tmp.val), &val, false);
//        if(AADC_UNLIKELY(recording)) CAAD_iVarAddAssign(&(other.val), &val);
#endif
        val = tmp.val;
        AADC_CHECK_VALUE(val);

        return *this;
    }
    idouble& operator /= (const idouble& other) noexcept {
        idouble tmp;
        // TODO: Implement *= more efficiently
        tmp.val = val / other.val;
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(recording)) CAAD_iVarDiv(&val, &other.val, &tmp.val);
        if(AADC_UNLIKELY(recording)) CAAD_iVarAssign(&(tmp.val), &val, false);
//        if(AADC_UNLIKELY(recording)) CAAD_iVarAddAssign(&(other.val), &val);
#endif
        val = tmp.val;
        AADC_CHECK_VALUE(val);

        return *this;
    }
    idouble operator -() const noexcept {
        return idouble(0.0) - (*this);
    }
    idouble operator +() const noexcept {
        return (*this);
    }
    aadc::AADCArgument markAsInput(int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsInput(&val, 0, i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsInput(const std::string& str_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsInput(&val, str_id.c_str(), i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCScalarArgument markAsScalarInput() const {
        CAAD_iVarMarkAsScalarInput(&val);
        return aadc::AADCScalarArgument(varIndex());
    }
    aadc::AADCArgument markAsDiff(int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsDiff(&val, 0, i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsDiff(const std::string& str_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsDiff(&val, str_id.c_str(), i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsInputNoDiff(int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsInputNoDiff(&val, 0, i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsInputNoDiff(const std::string& str_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsInputNoDiff(&val, str_id.c_str(), i1, i2, i3);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsArrayInput(uint32_t indx) const {
        CAAD_iVarMarkAsArrayInput(&val, indx);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument markAsArrayInputNoDiff(uint32_t indx) const {
        CAAD_iVarMarkAsArrayInputNoDiff(&val, indx);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCResult markAsOutput(int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsOutput(&val, 0, i1, i2, i3);
        return aadc::AADCResult(varIndex());
    }
    aadc::AADCResult markAsOutput(const std::string& str_id, int64_t i1 = -1, int64_t i2 = -1, int64_t i3 = -1) const {
        CAAD_iVarMarkAsOutput(&val, str_id.c_str(), i1, i2, i3);
        return aadc::AADCResult(varIndex());
    }
    aadc::AADCArgument forceAsVariable() const {
        CAAD_iVarForceVariable(&val);
        return aadc::AADCArgument(varIndex());
    }
    void forceAsDiff() const {
        CAAD_iVarForceDiff(&val);
    }
    aadc::AADCArgument initVarAtAddress(uint64_t addr) const {
        CAAD_iInitVarWithSpecificAddress(&val, false, addr);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument initScalarVarAtAddress(uint64_t addr) const {
        CAAD_iInitVarWithSpecificAddress(&val, true, addr);
        return aadc::AADCArgument(varIndex());
    }
    aadc::AADCArgument saveToStack() const {
        CAAD_iVarSaveToStack(&val);
        return aadc::AADCArgument(varIndex());
    }
    uint64_t varAddress() const noexcept {
        return CAAD_iVarIndex(&val);
    }
public:
    static void CheckPoint() {
        if(AADC_UNLIKELY(recording)) CAAD_CheckPoint();
    }
    static void DoNotZeroInitAllCurrentDiffVariables() {
        if(AADC_UNLIKELY(recording)) CAAD_DoNotZeroInitAllCurrentDiffVariables();
    }
    static void LoopPulse(const uint64_t i) {
        if(AADC_UNLIKELY(recording)) CAAD_LoopPulse(i);
    }
public: // for boost
    //! RUI
    operator double() const {return this->val;} 
    //! RUI
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & val;
    }
private:
    uint64_t varIndex() const noexcept {
        return CAAD_iVarIndex(&val);
    }
public:
    double val;
    //! COMENTEI ISTO
//operator double() const ;
//! COMENTEI ISTO
};


inline idouble operator + (const idouble& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val + b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarAdd(&a.val, &b.val, &c.val, false, false);
#endif
    return c;
}
inline idouble operator + (const idouble&& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val + b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarAdd(&a.val, &b.val, &c.val, true, false);
#endif
    return c;
}
inline idouble operator + (const idouble& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val + b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarAdd(&a.val, &b.val, &c.val, false, true);
#endif
    return c;
}
inline idouble operator + (const idouble&& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val + b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarAdd(&a.val, &b.val, &c.val, true, true);
#endif
    return c;
}

#define AADC_MIX_BINARY_OP(op, type, itype) \
inline idouble operator op (type a, itype b) noexcept { \
    return idouble(a) op b; \
} \
inline idouble operator op (itype a, type b) noexcept { \
    return a op idouble(b); \
}

#define AADC_MIX_BINARY_OP_ALL(op, itype) \
AADC_MIX_BINARY_OP(op, double, itype) \
AADC_MIX_BINARY_OP(op, int, itype) \
AADC_MIX_BINARY_OP(op, unsigned int, itype) \
AADC_MIX_BINARY_OP(op, int64_t, itype) \
AADC_MIX_BINARY_OP(op, uint64_t, itype)

AADC_MIX_BINARY_OP_ALL(+, const idouble&)
AADC_MIX_BINARY_OP_ALL(+, const idouble&&)

inline idouble operator - (const idouble& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val - b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarSub(&a.val, &b.val, &c.val, false, false);
#endif
    return c;
}
inline idouble operator - (const idouble&& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val - b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarSub(&a.val, &b.val, &c.val, true, false);
#endif
    return c;
}
inline idouble operator - (const idouble& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val - b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarSub(&a.val, &b.val, &c.val, false, true);
#endif
    return c;
}
inline idouble operator - (const idouble&& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val - b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarSub(&a.val, &b.val, &c.val, true, true);
#endif
    return c;
}

AADC_MIX_BINARY_OP_ALL(-, const idouble&)
AADC_MIX_BINARY_OP_ALL(-, const idouble&&)

inline idouble operator * (const idouble& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val * b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarMult(&a.val, &b.val, &c.val, false, false);
#endif
    return c;
}
inline idouble operator * (const idouble&& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val * b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarMult(&a.val, &b.val, &c.val, true, false);
#endif
    return c;
}
inline idouble operator * (const idouble& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val * b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarMult(&a.val, &b.val, &c.val, false, true);
#endif
    return c;
}
inline idouble operator * (const idouble&& a, const idouble&& b) noexcept {
    idouble c;
    c.val = a.val * b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarMult(&a.val, &b.val, &c.val, true, true);
#endif
    return c;
}

AADC_MIX_BINARY_OP_ALL(*, const idouble&)
AADC_MIX_BINARY_OP_ALL(*, const idouble&&)


inline idouble operator / (const idouble& a, const idouble& b) noexcept {
    idouble c;
    c.val = a.val / b.val;
    AADC_CHECK_VALUE(c.val);

#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarDiv(&a.val, &b.val, &c.val);
#endif
    return c;
}

AADC_MIX_BINARY_OP_ALL(/, const idouble&)



inline idouble& operator++(idouble& a) noexcept { a+=1.0; return a; }
inline idouble& operator--(idouble& a) noexcept { a-=1.0; return a; }
inline idouble operator++(idouble& a, int) noexcept { idouble olda(a); a+=1.0; return olda; }
inline idouble operator--(idouble& a, int) noexcept { idouble olda(a); a-=1.0; return olda; }

#define AADC_OPS_NAMESPACE std

#ifdef AADC_IDOUBLE_ACTIVE  
#define IDOUBLE_UNARY_FUNC(name, iname) \
extern "C" AADC_API void CAAD_iVar ## iname(const double* in, double* out); \
namespace AADC_OPS_NAMESPACE { \
inline idouble name(const idouble& a) noexcept { \
    idouble c; \
    c.val = name(a.val); \
    AADC_CHECK_VALUE(c.val); \
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVar ## iname(&a.val, &c.val); \
    return c; \
} \
} 
#else
#define IDOUBLE_UNARY_FUNC(name, iname) \
extern "C" void CAAD_iVar ## iname(const double* in, double* out); \
namespace AADC_OPS_NAMESPACE { \
inline idouble name(const idouble& a) noexcept { \
    idouble c; \
    c.val = name(a.val); \
    return c; \
} \
} 
#endif

// IDOUBLE_UNARY_FUNC(log, Log)
IDOUBLE_UNARY_FUNC(log10, Log10)
IDOUBLE_UNARY_FUNC(sin, Sin)
IDOUBLE_UNARY_FUNC(cos, Cos)
IDOUBLE_UNARY_FUNC(tan, Tan)
IDOUBLE_UNARY_FUNC(asin, Asin)
IDOUBLE_UNARY_FUNC(acos, Acos)
IDOUBLE_UNARY_FUNC(atan, Atan)
IDOUBLE_UNARY_FUNC(sinh, Sinh)
IDOUBLE_UNARY_FUNC(cosh, Cosh)

// IDOUBLE_UNARY_FUNC(exp, Exp)
IDOUBLE_UNARY_FUNC(sqrt, Sqrt)
IDOUBLE_UNARY_FUNC(tanh, Tanh)
IDOUBLE_UNARY_FUNC(ceil, Ceil)
IDOUBLE_UNARY_FUNC(floor, Floor)

#if _MSC_VER == 1700
static inline double log2(double n)
{
    return log(n) * FM_DOUBLE_LOG2OFE;
}

static inline double exp2(double n)
{
    return exp(n * (1.0 / FM_DOUBLE_LOG2OFE));
}
#endif

IDOUBLE_UNARY_FUNC(log2, Log2)
IDOUBLE_UNARY_FUNC(exp2, Exp2)
	
#if !_WIN64 || _MSC_VER > 1700
IDOUBLE_UNARY_FUNC(expm1, Expm1)
IDOUBLE_UNARY_FUNC(log1p, Log1p)
IDOUBLE_UNARY_FUNC(asinh, Asinh)
IDOUBLE_UNARY_FUNC(acosh, Acosh)
IDOUBLE_UNARY_FUNC(atanh, Atanh)
IDOUBLE_UNARY_FUNC(erf, Erf)
IDOUBLE_UNARY_FUNC(erfc, Erfc)
IDOUBLE_UNARY_FUNC(cbrt, Cbrt)
IDOUBLE_UNARY_FUNC(round, Round)
IDOUBLE_UNARY_FUNC(trunc, Trunc)
// IDOUBLE_UNARY_FUNC(rint, Rint)
// IDOUBLE_UNARY_FUNC(nearbyint, Nearbyint)
#endif

#ifdef AADC_IDOUBLE_ACTIVE  
#define IDOUBLE_BINARY_FUNC(name, iname) \
extern "C" AADC_API void CAAD_iVar ## iname(const double* in1, const double* in2, double* out); \
namespace AADC_OPS_NAMESPACE { \
inline idouble name(const idouble& a, const idouble& b) noexcept { \
    idouble c; \
    c.val = name(a.val, b.val); \
    AADC_CHECK_VALUE(c.val); \
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVar ## iname(&a.val, &b.val, &c.val); \
    return c; \
} \
} 
#else
#define IDOUBLE_BINARY_FUNC(name, iname) \
extern "C" void CAAD_iVar ## iname(const double* in1, const double* in2, double* out); \
namespace AADC_OPS_NAMESPACE { \
inline idouble name(const idouble& a, const idouble& b) noexcept { \
    idouble c; \
    c.val = name(a.val, b.val); \
    return c; \
} \
} 
#endif

IDOUBLE_BINARY_FUNC(pow, Pow)
IDOUBLE_BINARY_FUNC(atan2, Atan2)
IDOUBLE_BINARY_FUNC(fmod, Fmod)

namespace AADC_OPS_NAMESPACE {
inline idouble exp(const idouble& x) noexcept {
    return AADC_OPS_NAMESPACE::exp2(FM_DOUBLE_LOG2OFE * x);
}

inline idouble log(const idouble& x) noexcept {
    return AADC_OPS_NAMESPACE::log2(x) * (1.0 / FM_DOUBLE_LOG2OFE);
}

inline idouble pow(const idouble& x, double y) noexcept {
    return AADC_OPS_NAMESPACE::pow(x, idouble(y));
}
inline idouble pow(double x, const idouble& y) noexcept {
    return AADC_OPS_NAMESPACE::pow(idouble(x), y);
}

inline idouble fmod(const idouble& x, double y) noexcept {
    return AADC_OPS_NAMESPACE::fmod(x, idouble(y)); 
}
inline idouble fmod(double x, const idouble& y) noexcept {
    return AADC_OPS_NAMESPACE::fmod(idouble(x), y); 
}

inline double cdf_normal(double x) noexcept {
  static const double RT2PI = sqrt(4.0*acos(0.0));

  static const double SPLIT = 7.07106781186547;

  static const double N0 = 220.206867912376;
  static const double N1 = 221.213596169931;
  static const double N2 = 112.079291497871;
  static const double N3 = 33.912866078383;
  static const double N4 = 6.37396220353165;
  static const double N5 = 0.700383064443688;
  static const double N6 = 3.52624965998911e-02;
  static const double M0 = 440.413735824752;
  static const double M1 = 793.826512519948;
  static const double M2 = 637.333633378831;
  static const double M3 = 296.564248779674;
  static const double M4 = 86.7807322029461;
  static const double M5 = 16.064177579207;
  static const double M6 = 1.75566716318264;
  static const double M7 = 8.83883476483184e-02;

  const double z = fabs(x);
  double c = 0.0;

  if(z<=37.0)
  {
    const double e = exp(-z*z/2.0);
    if(z<SPLIT)
    {
      const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
      const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
      c = e*n/d;
    }
    else
    {
      const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
      c = e/(RT2PI*f);
    }
  }
  return x<=0.0 ? c : 1-c;
}

inline double cdf_normal(const double& x, double& pdf) noexcept 
{
    pdf = exp(-x*x * 0.5) / 2.506628274631;
    return cdf_normal(x);
}

inline idouble cdf_normal(const idouble& x, idouble& pdf) noexcept 
{
    idouble res;
    res.val = cdf_normal(x.val, pdf.val);
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarCDFNorm(&x.val, &res.val, &pdf.val);
    return res;
}

};

using AADC_OPS_NAMESPACE::log10;
using AADC_OPS_NAMESPACE::sin;
using AADC_OPS_NAMESPACE::cos;
using AADC_OPS_NAMESPACE::tan;
using AADC_OPS_NAMESPACE::asin;
using AADC_OPS_NAMESPACE::acos;
using AADC_OPS_NAMESPACE::atan;
using AADC_OPS_NAMESPACE::sinh;
using AADC_OPS_NAMESPACE::cosh;

using AADC_OPS_NAMESPACE::sqrt;
using AADC_OPS_NAMESPACE::tanh;
using AADC_OPS_NAMESPACE::ceil;
using AADC_OPS_NAMESPACE::floor;

using AADC_OPS_NAMESPACE::log2;

using AADC_OPS_NAMESPACE::exp;
using AADC_OPS_NAMESPACE::log;
using AADC_OPS_NAMESPACE::pow;
using AADC_OPS_NAMESPACE::fmod;


inline void CAAD_Check(const idouble& val, const char* msg) {
  //  CAAD_Assert(msg, &val.val);
}

inline void CAAD_Check(const double& val, const char* msg) {
//    CAAD_Assert(msg, &val.val);
}

namespace aadc {
class ExtVarIndex {
public:
    ExtVarIndex(const idouble& var, bool save_to_stack = true, bool force_diff = true)
        : index(idouble::recording ? var.forceAsVariable() : aadc::AADCArgument())
    {
        if (idouble::recording) {
            if (save_to_stack) var.saveToStack();
            if (force_diff) var.forceAsDiff();
        }
    }
    ExtVarIndex(const ExtVarIndex& other)
        : index(other.index)
    {}
    operator uint64_t () const {
        return index.getIndex();
    }
private:
    const aadc::AADCArgument index;
};
};

typedef double mdouble;

#ifndef AAD_NAMESPACE_BEGIN
#define AAD_NAMESPACE_BEGIN
#endif

#ifndef AAD_NAMESPACE_END
#define AAD_NAMESPACE_END
#endif

// Helper functions. Use with caution

inline double passive_value(const idouble& v) {
    return v.val;
}
inline double value(const idouble& v) {
    return v.val;
}

inline std::ostream& operator<<(std::ostream& os, const idouble& a)
{
#ifdef AADC_IDOUBLE_ACTIVE            
    if (idouble::recording) {
        bool is_random(CAAD_iVarIsRandom(&(a.val)));
        bool is_diff(CAAD_iVarIsDiff(&(a.val)));
        if (is_random || is_diff) {
            os << "[AAD"
            << (is_random ? "[rv] " : "")
            << (is_diff ? "[adj] " : "")
            << ":" << a.varAddress() << "," << a.val << "]";
        } else {
            os << a.val;
        }
    }
    else 
#endif    
    os << a.val;
    return os;
}
namespace aadc {
namespace detail {
template<typename T>        
struct ActiveToPassiveMap {
    typedef T PassiveType;
};
template<>        
struct ActiveToPassiveMap<idouble> {
    typedef double PassiveType;
};

};
template<class T>
inline typename detail::ActiveToPassiveMap<T>::PassiveType extractPassiveValue(const T& val, const char * file, uint64_t line) {
    return val;
}
template<>
inline double extractPassiveValue<idouble>(const idouble& val, const char * file, uint64_t line) {
#ifdef AADC_IDOUBLE_ACTIVE  
    if(AADC_UNLIKELY(idouble::recording)) return CAAD_iVarExtractPassive(&(val.val), file, line);
#endif
    return val.val;
}

template<class mdouble>
inline mdouble markVariableAsRandomInput(const double value) {
    return value;
}

template<>
inline idouble markVariableAsRandomInput<idouble>(const double value) {
    idouble var(value);
#ifdef AADC_IDOUBLE_ACTIVE  
    if(AADC_UNLIKELY(idouble::recording)) {
        CAAD_iVarMarkVariableAsRandomInput(&(var.val));
    }
#endif
    return var;
}
}

#define AAD_PASSIVE(ival) aadc::extractPassiveValue(ival, __FILE__, __LINE__) 


namespace std {
    inline string to_string( const idouble& v ) {
        return to_string(v.val);
    }
    inline bool isfinite(const idouble& v) {
        return isfinite(AAD_PASSIVE(v));
    }

  template<>
    struct numeric_limits<idouble> : public numeric_limits<double> {};

}; 

#endif // CAAD_IDOUBLE_H
