#pragma once
#include <cmath>
#include <algorithm>


inline double iIf(const bool cond, const double a, const double b) {
    return cond ? a : b;
}
inline bool iIf(const bool cond, const bool a, const bool b) {
    return (cond && a) || (!cond && b);
}
inline int iIf(const bool& cond, const int& a, const int& b) {
    return cond ? a : b;
}

inline void condAssign(double& a, const bool cond, const double b) {
    a = iIf(cond, b, a);
}

namespace aadc {
    // template versions of math functions. Needed to help conversions with function pointers
    template<class T>
    inline T abs(const T v) {
        return std::abs(v);
    }
    template<class T>
    inline T sqrt(const T v) {
        return std::sqrt(v);
    }
    
    template<class T>
    inline T log(const T v) {
        return std::log(v);
    }
    
    template<class T>
    inline T exp(const T v) {
        return std::exp(v);
    }
    template<class T>
    inline T cos(const T v) {
        return std::cos(v);
    }
    template<class T>
    inline T sin(const T v) {
        return std::sin(v);
    }    
    template <class T>
    inline T sign (const T& z)
    {
        return iIf(z == 0.0, 0.0 , iIf(z < 0.0, -1.0 , 1.0));
    }
    template <class T>
    inline T signbit(const T& z)
    {
    return iIf(z < 0.0, 1.0 , 0.0);
    }    

    template <class T>
    inline T max (const T& a, const T& b)
    {
        return std::max(a,b);
    }
    template <class T>
    inline T min (const T& a, const T& b)
    {
        return std::min(a,b);
    }

    template<class T>
    inline T cdf_normal(T x)
    {
        return std::erfc(-x/std::sqrt(2))/2;
    }

};

namespace aadcBoolOps {};