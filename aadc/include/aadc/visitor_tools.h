#pragma once
#include <iomanip>
#include <aadc/aadc.h>
#include <aadc/aadc_matrix.h>

namespace aadc {

////////////////////////////////////////////////////
//
//  VectorType<T>::VecType
//
//  Type traits to map from element type to corresponding vector type.
//
////////////////////////////////////////////////////

template<class T>
struct VectorType {
    typedef std::vector<T> VecType;
};

template<>
struct VectorType<__m256d> {
    typedef mmVector<__m256d> VecType;
};

template<>
struct VectorType<__m512d> {
    typedef mmVector<__m512d> VecType;
};
template<>
struct VectorType<AADCArgument> {
    typedef VectorArg VecType;
};

template<>
struct VectorType<AADCResult> {
    typedef VectorRes VecType;
};

////////////////////////////////////////////////////
//
//  InitializeVisitor() 
//
//  Visitor structure intendent to be applied with xVAResults, xVADiffResuls
//  Initialize object
//
////////////////////////////////////////////////////
template<typename numberType>
struct InitializeVisitor {

    void visit (numberType& val, const aadc::AADCArgument& arg) const {
        val = aadc::mmSetConst<numberType>(0.);
    }

    void visit (numberType& val, const aadc::AADCResult& arg) const {
        val =aadc::mmSetConst<numberType>(0.);
    }
    
    void visit (mmVector<numberType>& val, const aadc::VectorArg& arg) const {
        val = mmVector<numberType>(arg.size(), aadc::mmSetConst<numberType>(0.));
    }

    void visit (mmVector<numberType>& val, const aadc::VectorRes& arg) const {
        val = mmVector<numberType>(arg.size(), aadc::mmSetConst<numberType>(0.));
    }

    void visit (std::vector<double>& val, const aadc::VectorArg& arg) const {
        val = std::vector<double>(arg.size(), 0.0);
    }

    void visit (std::vector<double>& val, const aadc::VectorRes& arg) const {
        val = std::vector<double>(arg.size(), 0.0);
    }

};

////////////////////////////////////////////////////
//
//  MMAccumulateVisitorWS
//
//  Visitor structure accumulates avx vector results taken from AADCWorkSpace
//  Implements <mmtype> +=
//
//  _ws     AADC work space
//
////////////////////////////////////////////////////

template<typename mmType>
class MMAccumulateVisitorWS {
public:
    MMAccumulateVisitorWS(const aadc::AADCWorkSpace<mmType>& _ws) : ws(_ws) {} 

    void visit(mmType& val, const aadc::AADCArgument& arg) const {
        val = aadc::mmAdd(val, ws.diff(arg));
    }
    void visit(mmType& val, const aadc::AADCResult& arg) const {
        val = aadc::mmAdd(val, ws.val(arg));
    }
    void visit(mmVector<mmType>& val, const aadc::VectorArg& arg) const {
        for (int i=0; i<arg.size(); i++) val[i] = aadc::mmAdd(val[i], ws.diff(arg[i]));
    }    
    void visit(mmVector<mmType>& val, const aadc::VectorRes& arg) const {
        for (int i=0; i<arg.size(); i++) val[i] = aadc::mmAdd(val[i], ws.val(arg[i]));
    }
private:
    const aadc::AADCWorkSpace<mmType>& ws;
};

////////////////////////////////////////////////////
//
//  MMSumReduceVisitor
// 
//  Implements a coordinateWise sum
//
//  _norm    normalization coefficient 
//
////////////////////////////////////////////////////

template<typename mmType>
class MMSumReduceVisitor {
public:
    MMSumReduceVisitor(const double& _norm) : norm(_norm) {}

    void visit(double& val, const mmType& mm_val) const {
        val = aadc::mmSum(mm_val)*norm;
    }
    void visit(std::vector<double>& vec, const mmVector<mmType>& mm_vec) const {
        for (int i=0; i<vec.size(); i++) vec[i] = aadc::mmSum(mm_vec[i])*norm;
    }
private:
    double norm;
};  

////////////////////////////////////////////////////
//
//  AccumulateVisitor
//
//  Collects information from various threads 
//
////////////////////////////////////////////////////

template<typename T>
class AccumulateVisitor {
public:
    void visit(double& val, const double& val2) const {
        val += val2;
    }
    void visit(std::vector<double>& vec, const std::vector<double>& vec2) const {
        std::transform (vec.begin(), vec.end(), vec2.begin(), vec.begin(), std::plus<double>());
    }
};  

}; // namespace aadc