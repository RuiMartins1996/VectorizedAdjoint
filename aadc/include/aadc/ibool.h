#ifndef CAAD_IBOOL_H
#define CAAD_IBOOL_H
#include <csignal>

#include "aadc/idouble.h"

extern "C" AADC_API void CAAD_iBoolVarConstructor(uint64_t*);
extern "C" AADC_API void CAAD_iBoolVarConstructorConstant(uint64_t*, const bool*);

extern "C" AADC_API void CAAD_iBoolVarDestructor(uint64_t*);

extern "C" AADC_API void CAAD_iBoolVarAssign(const uint64_t* in, uint64_t* out);

extern "C" AADC_API void CAAD_iVarAnd(const uint64_t* in1, const uint64_t* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarOr(const uint64_t* in1, const uint64_t* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarXor(const uint64_t* in1, const uint64_t* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarNot(const uint64_t* in, uint64_t* out);

extern "C" AADC_API void CAAD_iVarLess(const double* in1, const double* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarLessEqual(const double* in1, const double* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarNotEqual(const double* in1, const double* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarEqual(const double* in1, const double* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarGreaterEqual(const double* in1, const double* in2, uint64_t* out);
extern "C" AADC_API void CAAD_iVarGreater(const double* in1, const double* in2, uint64_t* out);

extern "C" AADC_API uint64_t CAAD_iVarExtractPassiveBool(const uint64_t* var, const char* loc_msg, uint64_t line);

extern "C" AADC_API void CAAD_iVarIf(const uint64_t* cond, const double* in1, const double* in2, double* out);

extern "C" AADC_API bool CAAD_iBoolVarIsRandom(const uint64_t* v);

class ibool {
public:
    ibool(const bool& c) : val(c) { 
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(idouble::recording)) CAAD_iBoolVarConstructorConstant(&val, &c);
#endif
    }
    ibool() {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(idouble::recording)) CAAD_iBoolVarConstructor(&val);
#endif
    }
    ibool(const ibool& other) : val(other.val) {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(idouble::recording)) CAAD_iBoolVarAssign(&other.val, &val);
#endif
    }
    ~ibool() {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(idouble::recording)) CAAD_iBoolVarDestructor(&val);
#endif
    }

    ibool& operator = (const ibool& other) {
#ifdef AADC_IDOUBLE_ACTIVE        
        if(AADC_UNLIKELY(idouble::recording)) CAAD_iBoolVarAssign(&(other.val), &val);
#endif
        val = other.val;
        return *this;
    }
#ifdef AADC_ALLOW_TO_PASSIVE_BOOL
    operator bool() const { 
#ifdef AADC_IDOUBLE_ACTIVE        
        if(idouble::recording && CAAD_iBoolVarIsRandom(&val)) {
            CAAD_iVarExtractPassiveBool(&val, "ibool -> bool at UNKNOWN Location. Use AADC_BreakOnActiveBoolConversion option", 0);
        }
#endif
        return val;
    }
#endif
    // TODO: Remove from production code
    uint64_t varIndex() {
        return CAAD_iVarIndex((double*)&val);
    }

public:
    uint64_t val;
};
namespace aadcBoolOps {
inline ibool operator && (const ibool& a, const ibool& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarAnd(&a.val, &b.val, &c.val);
#endif
    c.val = a.val && b.val;
    return c;
}
inline ibool operator && (bool a, const ibool& b) {
    return ibool(a) && b;
}
inline ibool operator && (const ibool& a, const bool& b) {
    return a && ibool(b);
}

inline ibool operator || (const ibool& a, const ibool& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarOr(&a.val, &b.val, &c.val);
#endif
    c.val = a.val || b.val;
    return c;
}
inline ibool operator || (const bool& a, const ibool& b) {
    return ibool(a) || b;
}
inline ibool operator || (const ibool& a, const bool& b) {
    return ibool(b) || a;
}
};

inline ibool operator != (const ibool& a, const ibool& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarXor(&a.val, &b.val, &c.val);
#endif
    c.val = a.val != b.val;
    return c;
}
inline ibool operator != (const bool& a, const ibool& b) {
    return ibool(a) != b;
}
inline ibool operator != (const ibool& a, const bool& b) {
    return (a) != ibool(b);
}

inline ibool operator!(const ibool& a) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarNot(&a.val,  &c.val);
#endif
    c.val = !a.val;                                         
    return c;
}
inline ibool operator < (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarLess(&a.val, &b.val, &c.val);
#endif
    c.val = a.val < b.val;
    return c;
}

#define AADC_MIX_CMP_BINARY_OP(op, type) \
inline ibool operator op (type a, const idouble& b) { \
    return idouble(a) op b; \
} \
inline ibool operator op (const idouble& a, type b) { \
    return a op idouble(b); \
}

#define AADC_MIX_CMP_BINARY_OP_ALL(op) \
AADC_MIX_CMP_BINARY_OP(op, double) \
AADC_MIX_CMP_BINARY_OP(op, int) \
AADC_MIX_CMP_BINARY_OP(op, unsigned int) \
AADC_MIX_CMP_BINARY_OP(op, int64_t) \
AADC_MIX_CMP_BINARY_OP(op, uint64_t)

AADC_MIX_CMP_BINARY_OP_ALL(<)

inline ibool operator <= (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarLessEqual(&a.val, &b.val, &c.val);
#endif
    c.val = a.val <= b.val;
    return c;
}

AADC_MIX_CMP_BINARY_OP_ALL(<=)

inline ibool operator != (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarNotEqual(&a.val, &b.val, &c.val);
#endif
    c.val = a.val != b.val;
    return c;
}

AADC_MIX_CMP_BINARY_OP_ALL(!=)

inline ibool operator == (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarEqual(&a.val, &b.val, &c.val);
#endif
    c.val = a.val == b.val;
    return c;
}

AADC_MIX_CMP_BINARY_OP_ALL(==)

inline ibool operator >= (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarGreaterEqual(&a.val, &b.val, &c.val);
#endif
    c.val = a.val >= b.val;
    return c;
}

AADC_MIX_CMP_BINARY_OP_ALL(>=)

inline ibool operator > (const idouble& a, const idouble& b) {
    ibool c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarGreater(&a.val, &b.val, &c.val);
#endif
    c.val = a.val > b.val;
    return c;
}

AADC_MIX_CMP_BINARY_OP_ALL(>)

inline idouble iIf(const ibool& cond, const idouble& a, const idouble& b) {
    idouble c;
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) CAAD_iVarIf(&cond.val, &a.val, &b.val, &c.val);
#endif
    c.val = cond.val ? a.val : b.val;
    return c;
}

inline idouble iIf(const ibool& cond, const double a, const double b) {
    return iIf(cond, idouble(a), idouble(b));
}

inline idouble iIf(const ibool& cond, const idouble& a, const double b) {
    return iIf(cond, (a), idouble(b));
}

inline idouble iIf(const ibool& cond, const double a, const idouble& b) {
    return iIf(cond, idouble(a), (b));
}

inline ibool iIf(const ibool& cond, const ibool& a, const ibool& b) {
    using namespace aadcBoolOps;
    return (cond && a) || (!cond && b);
}
inline ibool iIf(const ibool& cond, const bool a, const ibool& b) {
    using namespace aadcBoolOps;
    return (cond && a) || (!cond && b);
}
inline ibool iIf(const ibool& cond, const ibool& a, const bool b) {
    using namespace aadcBoolOps;
    return (cond && a) || (!cond && b);
}

inline void condAssign(idouble& a, const ibool& cond, const idouble& b) {
    a = iIf(cond, b, a);
}

namespace aadc {
namespace detail {
template<>        
struct ActiveToPassiveMap<ibool> {
    typedef bool PassiveType;
};

};

template<>
inline bool extractPassiveValue<ibool>(const ibool& val, const char * file, uint64_t line) {
#ifdef AADC_IDOUBLE_ACTIVE        
    if(AADC_UNLIKELY(idouble::recording)) return CAAD_iVarExtractPassiveBool(&(val.val), file, line);
#endif
    return val.val;
}
}


namespace std {
    inline idouble max(const idouble& a,const idouble& b) {
        return iIf(a < b, b, a);
    }
    inline idouble max(const idouble& a,const double b) {
        return iIf(a < b, b, a);
    }
    inline idouble max(const double a,const idouble& b) {
        return iIf(a < b, b, a);
    }
    inline idouble min(const idouble& a,const idouble& b) {
        return iIf(a < b, a, b);
    }
    inline idouble min(const idouble& a,const double b) {
        return iIf(a < b, a, b);
    }
    inline idouble min(const double a,const idouble& b) {
        return iIf(a < b, a, b);
    }
};

// Common ops

//inline idouble abs (const idouble& a) {
//    return std::max(a, -a);
//}

namespace std {
    inline idouble abs(const idouble& x) {
        return std::max(x, -x);
    }
    inline idouble fabs(const idouble& x) {
        return std::max(x, -x);
    }

    inline idouble copysign(const idouble& a, const idouble& b) {
        return iIf(b>=0.0, a, -a);
    }
};

using std::fabs;

#include <aadc/aadc_compat.h>

#endif // CAAD_IBOOL_H
