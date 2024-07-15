#pragma once
#include <aadc/idouble.h>
#include <aadc/idouble.h>
#include <Eigen/Dense>

namespace Eigen {
template<> struct NumTraits<idouble>
    : NumTraits<double> 
{
    typedef idouble Real;
    typedef idouble NonInteger;
    typedef idouble Nested;
    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3
    };
};

template<> struct NumTraits<::std::complex<idouble>>
 : NumTraits<::std::complex<double>> 
{
typedef idouble Real;
typedef ::std::complex<idouble> NonInteger;
typedef ::std::complex<idouble> Nested;
enum {
    IsComplex = 1,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
};
};

template<typename BinaryOp> struct ScalarBinaryOpTraits<double, idouble, BinaryOp>
{
  typedef idouble ReturnType;
};

template<typename BinaryOp> struct ScalarBinaryOpTraits<idouble, double, BinaryOp>
{
  typedef idouble ReturnType;
};
}
