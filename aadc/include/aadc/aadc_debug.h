#pragma once
#include "aadc/aadc_ext.h"
#include "aadc/idouble.h"

#include <iostream>
#include <sstream>
#include <iomanip>

extern "C" AADC_API void AADDebugStart(double adjoint_print_scaler = 1.0, int aadc_debug_avx_index_ = 0);

extern "C" AADC_API bool AADDebugNeedBumpBase();
extern "C" AADC_API void AADDebugStop();
extern "C" AADC_API void AADDebugResult(double& res);

extern "C" AADC_API double AADCDebugAdj(const double& x, const char *var_name, const char* file_name, int line_num, double bump_h);
idouble AADDebugResult(idouble res);

extern "C" AADC_API void CAAD_CPUState(
    bool forward, 
    const char* name0 = 0, const double* var0 = 0,
    const char* name1 = 0, const double* var1 = 0,
    const char* name2 = 0, const double* var2 = 0,
    const char* name3 = 0, const double* var3 = 0,
    const char* name4 = 0, const double* var4 = 0,
    const char* name5 = 0, const double* var5 = 0,
    const char* name6 = 0, const double* var6 = 0,
    const char* name7 = 0, const double* var7 = 0,
    const char* name8 = 0, const double* var8 = 0,
    const char* name9 = 0, const double* var9 = 0
);

extern "C" AADC_API void AADCDebugInit();

extern "C" AADC_API void CAAD_BinaryBreakPoint();

extern "C" AADC_API void AADCDebugPutRandomFlag(uint64_t counter, bool is_random);
extern "C" AADC_API bool AADCDebugGetRandomFlag(uint64_t counter);
extern "C" AADC_API void AADCDebugPutDiffFlag(uint64_t counter, bool is_diff);
extern "C" AADC_API bool AADCDebugGetDiffFlag(uint64_t counter);

extern "C" AADC_API uint64_t* AADCGetDebugVarCounter();
extern "C" AADC_API void AADCGetDebugVarCounterInc();
extern "C" AADC_API uint64_t AADCGetTotalDebugPrints();
extern "C" AADC_API uint64_t* AADCGetDebugPrintAVXIndex();
extern "C" AADC_API uint64_t AADCGetStopRecordingAt();
extern "C" AADC_API void AADCSetStopRecordingAt(uint64_t counter);

extern "C" AADC_API void AADCDebugFwdBumpVar(uint64_t var_index, double h);
extern "C" AADC_API void AADCDebugGetFwdBumpVar(uint64_t& var_index, double& h);
extern "C" AADC_API bool AADCDebugIsBumpFwdRun();

extern "C" AADC_API void AADCDebugPutRecValue(uint64_t counter, double val);
extern "C" AADC_API double AADCDebugGetRecValue(uint64_t counter);
extern "C" AADC_API void AADCDebugPutFwdValue(uint64_t counter, double val);
extern "C" AADC_API double AADCDebugGetFwdValue(uint64_t counter);
extern "C" AADC_API void AADCDebugPutAdjValue(uint64_t counter, double val);
extern "C" AADC_API double AADCDebugGetAdjValue(uint64_t counter);

extern "C" AADC_API void AADCDebugPutLineNum(uint64_t counter, int val);
extern "C" AADC_API int AADCDebugGetLineNum(uint64_t counter);

extern "C" AADC_API void AADCDebugPutFileName(uint64_t counter, const char* val);
extern "C" AADC_API const char* AADCDebugGetFileName(uint64_t counter);

extern "C" AADC_API bool AADCIsDebugPrintNow();
extern "C" AADC_API double AADCDebugPrintScaler();
extern "C" AADC_API const char* AADCDebugVarAttrString(const bool is_random, const bool is_adjoint );

extern "C" AADC_API void AADCDebugPutVarName(uint64_t counter, const char* val);
extern "C" AADC_API const char* AADCDebugGetVarName(uint64_t counter);

extern "C" AADC_API const char* AADCDebugGetCurrentBlockName();
extern "C" AADC_API void AADCDebugPutCurrentBlockName(const char * block_name);
extern "C" AADC_API const char* AADCDebugGetCurrentFileName();
extern "C" AADC_API void AADCDebugPutCurrentFileName(const char * file_name);
extern "C" AADC_API int AADCDebugGetCurrentFileLine();
extern "C" AADC_API void AADCDebugPutCurrentFileLine(int file_line);

extern "C" AADC_API void AADCDebugPrintOn();
extern "C" AADC_API void AADCDebugPrintOff();

namespace aadc {

class DebugPrintExtFuncWrapper : public aadc::ConstStateExtFunc {
public:
  DebugPrintExtFuncWrapper(double* _res_val, const double* _x_val, const char *var_name_, const char* file_name_, int line_num_)
    : is_random(CAAD_iVarIsRandom(_x_val))
    , is_diff(CAAD_iVarIsDiff(_x_val))
    , xi((is_random || is_diff) ? CAAD_iVarIndex(_x_val) : 0)
    , x_val(*_x_val)
    , counter(*AADCGetDebugVarCounter())
    , var_name(var_name_)
    , file_name(file_name_)
    , line_num(line_num_)
  {
    if (counter == AADCGetStopRecordingAt()) {
      #ifdef _MSC_VER
        __debugbreak();
      #else
      #include <signal.h>
        raise(SIGTRAP);
      #endif
    }
    *_res_val = *_x_val;
    if (is_random) {
      CAAD_iVarForceVariable(_x_val);
      CAAD_iVarForceVariable(_res_val);
    }
    if (is_diff) {
      CAAD_iVarForceDiff(_res_val);
      CAAD_iVarForceDiff(_x_val);
    }
    if (is_random || is_diff) {
      resi = CAAD_iVarIndex(_res_val);
    }

    if (AADCIsDebugPrintNow() || !std::isfinite(*_x_val)) {
      std::stringstream sstr;
      sstr << std::setprecision(15)
          << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') << counter 
          << "::Val:aadc_rec\t" 
          << *_x_val << " "
          << "v[" << xi << "]"
          << "{" << var_name << "}"
          << AADCDebugVarAttrString(is_random, is_diff)
          << "@ " << file_name << ":" << line_num
          << std::endl
      ;
      std::cout << sstr.str();
    }
    AADCDebugPutRecValue(counter, x_val);
    AADCDebugPutRandomFlag(counter, is_random);
    AADCDebugPutDiffFlag(counter, is_diff);
    AADCDebugPutLineNum(counter, line_num_);
    AADCDebugPutFileName(counter, file_name_);
    AADCDebugPutVarName(counter, var_name_);

    AADCGetDebugVarCounterInc();
  }

  template<typename mmType>
  void forward(mmType* v) const {
    if (counter == AADCGetStopRecordingAt()) {
      #ifdef _MSC_VER
        __debugbreak();
      #else
      #include <signal.h>
        raise(SIGTRAP);
      #endif
    }
    if (!is_random) {
      if (AADCIsDebugPrintNow()) {
      std::stringstream sstr;
      sstr << std::setprecision(15)
          << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') << counter 
          << "::Val:aadc_val\t" 
          << x_val << " "
          << "{" << var_name << "}"
          << AADCDebugVarAttrString(is_random, is_diff)
          << " @ " << file_name << ":" << line_num
          << std::endl;
          std::cout << sstr.str();
      }
      if (!AADCDebugIsBumpFwdRun()) AADCDebugPutFwdValue(counter, x_val);
      return ;
    }
    v[resi] = v[xi];
    uint64_t bump_var_index;
    double bump_h;

    AADCDebugGetFwdBumpVar(bump_var_index, bump_h);
    if (bump_var_index == counter) {
      aadc::toDblPtr(v[resi])[*AADCGetDebugPrintAVXIndex()] += bump_h;
    }
    bool is_bumped_run(AADCDebugIsBumpFwdRun());
    if (!is_bumped_run) {
      AADCDebugPutFwdValue(counter, aadc::toDblPtr(v[xi])[*AADCGetDebugPrintAVXIndex()]);
    }
    if (AADCIsDebugPrintNow()) {
     std::stringstream sstr;
     sstr << std::setprecision(15)
        << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') << counter 
        << "::Val:aadc_val\t" 
        << toDblPtr(v[resi])[*AADCGetDebugPrintAVXIndex()] << " "
        << "v[" << xi << "]"
        << "{" << var_name << "}"
        << AADCDebugVarAttrString(is_random, is_diff)
        << " @ " << file_name << ":" << line_num
        << std::endl;
        std::cout << sstr.str();
    }
  }
  template<class mmType>
  void reverse(const mmType *v, mmType *d) const {
    if (!is_diff) {
      if (AADCIsDebugPrintNow()) {
        std::stringstream sstr;
        sstr << std::setprecision(15)
            << "AADCDebug::" << std::right <<std::setw(10) << std::setfill('0') << counter 
            << "::Adj:aadc_adj\t" 
            << "???"
            << " {" << var_name << "}" 
            << AADCDebugVarAttrString(is_random, is_diff)
            << "@ " << file_name << ":" << line_num
            << std::endl
        ;
        std::cout << sstr.str();
      }
      AADCDebugPutAdjValue(counter,0.);
      return ;
    }
    d[xi] = aadc::mmAdd(d[xi], d[resi]);
    if (AADCIsDebugPrintNow()) {
      std::stringstream sstr;
      sstr << std::setprecision(15)
          << "AADCDebug::" << std::right <<std::setw(10) << std::setfill('0') << counter 
          << "::Adj:aadc_adj\t" 
          << AADCDebugPrintScaler() * toDblPtr(d[resi])[*AADCGetDebugPrintAVXIndex()]
          << " {" << var_name << "}" 
          << AADCDebugVarAttrString(is_random, is_diff)
          << "@ " << file_name << ":" << line_num
          << std::endl
      ;
      std::cout << sstr.str();
    }
    AADCDebugPutAdjValue(counter, AADCDebugPrintScaler() * toDblPtr(d[resi])[*AADCGetDebugPrintAVXIndex()]);
  }

private:
  const bool is_random;
  const bool is_diff;
  const uint64_t xi;
  const double x_val;
  uint64_t resi;

  const int counter;
  const char *var_name;
  const char* file_name;
  const int line_num;
};

};

inline idouble AADCDebugAdj(const idouble& x, const char *var_name, const char* file_name, int line_num, double bump_h) {
    if (!idouble::recording) 
    {
      return AADCDebugAdj(x.val, var_name, file_name, line_num, bump_h);
    }
    idouble res;

    aadc::addConstStateExtFunction(
        std::make_shared<aadc::DebugPrintExtFuncWrapper>(
            &res.val, &x.val, var_name, file_name, line_num
        )
    );

    return res;
}

#define AADC_PRINT(var_name) \
AADCDebugAdj(var_name, #var_name, __FILE__, __LINE__, 0.00001)



inline void AADC_CPUState(bool forward, 
    const char* name0 = 0, const double* var0 = 0,
    const char* name1 = 0, const double* var1 = 0,
    const char* name2 = 0, const double* var2 = 0,
    const char* name3 = 0, const double* var3 = 0,
    const char* name4 = 0, const double* var4 = 0,
    const char* name5 = 0, const double* var5 = 0,
    const char* name6 = 0, const double* var6 = 0,
    const char* name7 = 0, const double* var7 = 0,
    const char* name8 = 0, const double* var8 = 0,
    const char* name9 = 0, const double* var9 = 0
) {

}

inline void AADC_CPUState(bool forward, 
    const char* name0 = 0, const idouble* var0 = 0,
    const char* name1 = 0, const idouble* var1 = 0,
    const char* name2 = 0, const idouble* var2 = 0,
    const char* name3 = 0, const idouble* var3 = 0,
    const char* name4 = 0, const idouble* var4 = 0,
    const char* name5 = 0, const idouble* var5 = 0,
    const char* name6 = 0, const idouble* var6 = 0,
    const char* name7 = 0, const idouble* var7 = 0,
    const char* name8 = 0, const idouble* var8 = 0,
    const char* name9 = 0, const idouble* var9 = 0
) {
  CAAD_CPUState(forward,
    name0, (double*)var0,
    name1, (double*)var1,
    name2, (double*)var2,
    name3, (double*)var3,
    name4, (double*)var4,
    name5, (double*)var5,
    name6, (double*)var6,
    name7, (double*)var7,
    name8, (double*)var8,
    name9, (double*)var9    
  );
}

#define AADC_CPU_STATE(forward, var1, var2, var3) \
AADC_CPUState(forward, \
#var1, &(var1), \
#var2, &(var2), \
#var3, &(var3) \
)

namespace aadc {


template<class mmType>
double computeFDDeriv(
  mmType central_value, aadc::AADCWorkSpace<mmType>& ws, const aadc::AADCFunctions<mmType>& aadc_func,
  int64_t var_index, aadc::AADCResult result,
  double h, const double* shifts, const double* weights, double div_by, int num_shifts
) {
  using namespace aadc;
  mmType sum(mmZero<mmType>());
  for (int i = 0; i < num_shifts; ++i) {
    AADCDebugFwdBumpVar(var_index, h * shifts[i]);
    aadc_func.forward(ws);
    mmType bumped_val = ws.val(result);

    sum = mmAdd(sum, mmMul(mmSetConst<mmType>( weights[i]) , mmSub(bumped_val,central_value)));
  }

  mmType deriv = mmDiv(sum,mmSetConst<mmType>(div_by * h));
  return toDblPtr(deriv)[*AADCGetDebugPrintAVXIndex()];
}


template<class stream, class mmType>
void debugTestAdjointsUsingFD(
  stream& logs
  , aadc::AADCWorkSpace<mmType>& ws
  , const aadc::AADCFunctions<mmType>& aadc_func
  , aadc::AADCResult result
  , bool check_fwd_vs_rec = false
  , bool print_all = false
) {
    AADCDebugPrintOff(); // disable printing during kernel execution for bump-and-revalue
  double central_shift[] = {-3.0,-2.0,   -1.0,  1.0,  2.0, 3.0};
  double central_w[] =     {-1.0, 9.0,  -45.0, 45.0, -9.0, 1.0};
  double central_div = 60.0;

  double left_shift[] = {-1.0,   -2.0,  -3.0,   -4.0,  -5.0,   -6.0};
  double left_w[] =     {72.0, -225.0, 400.0, -450.0, 360.0, -147.0};
  double left_div = 60.0;

  double right_shift[] = {  1.0,   2.0,   3.0,   4.0,    5.0,   6.0};
  double right_w[] =    { -72.0, 225.0,-400.0, 450.0, -360.0, 147.0};
  double right_div = 60.0;

  uint64_t total_print_vars(AADCGetTotalDebugPrints());

  mmType central_val = ws.val(result);
  int avx_i(*AADCGetDebugPrintAVXIndex());

  double eps(1e-5);
  
  // Check all forward values for sanity
  for (uint64_t print_i = 0; print_i < total_print_vars; ++print_i) {
    double var_val(AADCDebugGetFwdValue(print_i));
    bool is_random(AADCDebugGetRandomFlag(print_i));
    bool is_diff(AADCDebugGetDiffFlag(print_i));
    double var_rec_val(AADCDebugGetRecValue(print_i));
    if (!std::isfinite(var_val) || (check_fwd_vs_rec && (std::fabs(var_rec_val - var_val) > 1e-10))) {
      std::stringstream sstr;
      sstr << std::setprecision(15) << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Val[" << avx_i << "]:aadc_val\t" << var_val 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
//      logs << sstr.str();
      sstr << std::setprecision(15) << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Val[" << avx_i << "]:rec_val\t" << var_rec_val 
          << AADCDebugVarAttrString(is_random, is_diff)
      ;
      if (AADCDebugGetLineNum(print_i)) {
        sstr << "{" << AADCDebugGetVarName(print_i) << "} " << AADCDebugGetFileName(print_i) << ":" << AADCDebugGetLineNum(print_i);
      }
      sstr << std::endl;
      logs << sstr.str();
    }
  }
  for (int64_t print_i = int64_t(total_print_vars) - 1; print_i >= 0; --print_i) {
    bool is_random(AADCDebugGetRandomFlag(print_i));
    bool is_diff(AADCDebugGetDiffFlag(print_i));

    if (!is_random || !is_diff) {
      if (print_all) {
        std::stringstream sstr;

        sstr << std::setprecision(15) 
            << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
            << print_i << "::Adj:FD check N/A "
            << AADCDebugVarAttrString(is_random, is_diff);
        if (AADCDebugGetLineNum(print_i)) {
          sstr << " {" << AADCDebugGetVarName(print_i) << "} " << AADCDebugGetFileName(print_i) << ":" << AADCDebugGetLineNum(print_i);
        }

        sstr << std::endl;
        logs << sstr.str();
      }
      continue;
    }

    double var_val(AADCDebugGetFwdValue(print_i));

    double aadc_deriv = AADCDebugGetAdjValue(print_i);

    double h = std::max(fabs(var_val) * 0.0002,2e-5);

    double tol(std::max(fabs(aadc_deriv)*eps, 1e-5));

    double central_deriv = computeFDDeriv(central_val, ws, aadc_func, print_i, result, h, central_shift, central_w, central_div, sizeof(central_shift)/sizeof(central_shift[0]));
    bool central_deriv_ok = std::isfinite(central_deriv) && fabs(central_deriv - aadc_deriv) < tol;
    double left_deriv = computeFDDeriv(central_val, ws, aadc_func, print_i, result, h, left_shift, left_w, left_div, sizeof(left_shift)/sizeof(left_shift[0]));
    bool left_deriv_ok = std::isfinite(left_deriv) && fabs(left_deriv - aadc_deriv) < tol;
    double right_deriv = computeFDDeriv(central_val, ws, aadc_func, print_i, result, h, right_shift, right_w, right_div, sizeof(right_shift)/sizeof(right_shift[0]));
    bool right_deriv_ok = std::isfinite(right_deriv) && fabs(right_deriv - aadc_deriv) < tol;

    bool print_this = print_all ||
      !std::isfinite(var_val) || !std::isfinite(aadc_deriv) ||
      !(central_deriv_ok || left_deriv_ok || right_deriv_ok)
    ;

    if (print_this) {
      std::stringstream sstr;

      sstr << std::setprecision(15)      
          << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Adj[" << avx_i << "]:FD_mid \t" << central_deriv 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
      sstr << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Adj[" << avx_i << "]:FD_left\t" << left_deriv 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
      sstr << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Adj[" << avx_i << "]:FD_right\t" << right_deriv 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
      sstr << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Adj[" << avx_i << "]:aadc_adj\t" << aadc_deriv 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
      sstr << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Val[" << avx_i << "]:aadc_val\t" << var_val 
          << AADCDebugVarAttrString(is_random, is_diff)
          << std::endl
      ;
      sstr << "AADCDebug::" << std::right << std::setw(10) << std::setfill('0') 
          << print_i << "::Val[" << avx_i << "]:rec_val\t" << AADCDebugGetRecValue(print_i)
          << AADCDebugVarAttrString(is_random, is_diff)
      ;
        if (AADCDebugGetLineNum(print_i)) {
          sstr << "{" << AADCDebugGetVarName(print_i) << "} " << AADCDebugGetFileName(print_i) << ":" << AADCDebugGetLineNum(print_i);
        }
        sstr << std::endl
      ;
      
      logs << sstr.str();
    }
  }
}

class AADCDebugNamedCodeScope {
public:
  AADCDebugNamedCodeScope(
    const char * block_name,
    const char * file_name,
    int file_line
  ) 
    : prev_block_name(AADCDebugGetCurrentBlockName())
    , prev_file_name(AADCDebugGetCurrentFileName())
    , prev_file_line(AADCDebugGetCurrentFileLine())
  {
    AADCDebugPutCurrentBlockName(block_name);
    AADCDebugPutCurrentFileName(file_name);
    AADCDebugPutCurrentFileLine(file_line);
  }
  ~AADCDebugNamedCodeScope() {
    AADCDebugPutCurrentBlockName(prev_block_name);
    AADCDebugPutCurrentFileName(prev_file_name);
    AADCDebugPutCurrentFileLine(prev_file_line);
  }
private:
  const char* prev_block_name;
  const char* prev_file_name;
  int prev_file_line;
};


#define AADC_PRINT_DEBUG_COMMENT(block_name) \
AADCDebugNamedCodeScope aadc_debug_scope(block_name, __FILE__, __LINE__)


};