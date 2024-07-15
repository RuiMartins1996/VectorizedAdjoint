#ifndef AADC_EXT_H
#define AADC_EXT_H

#include <immintrin.h>
#include <vector>
#include <memory>

extern "C" AADC_API void CAAD_ExtFunction(const void* this_ptr, const void* fwd_ptr, const void* rev_ptr, const void* fwd_ptr_512, const void* rev_ptr_512);

namespace aadc {

class ConstStateExtFunc {
public:

// Derived class should implement the following two methods
//  template<typename mmType> void forward(mmType* v);
//  template<typename mmType> void reverse(const mmType* v, mmType *d);

  virtual ~ConstStateExtFunc() {}
};

inline std::vector<std::shared_ptr<ConstStateExtFunc> >& getRecordingConstStateExtFunctionRegistry() {
    static std::vector<std::shared_ptr<ConstStateExtFunc> > ext_funcs;
    return ext_funcs;
}

template<class FuncWrap, class mmType>
inline void callFuncWrapFwd(const FuncWrap* obj, mmType* v) {
  obj->template forward<mmType>(v);
}

template<class FuncWrap, class mmType>
inline void callFuncWrapRev(const FuncWrap* obj, mmType* v, mmType* d) {
  obj->template reverse<mmType>(v,d);
}

template<class FuncWrap>
void addConstStateExtFunction(const std::shared_ptr<FuncWrap>& func_wrap) {
  if (!idouble::recording) return ;
  getRecordingConstStateExtFunctionRegistry().push_back(func_wrap);

    auto obj_ptr = func_wrap.get();
    auto call_256_fwd_ptr = &callFuncWrapFwd<FuncWrap, __m256d>;
    auto call_256_rev_ptr = &callFuncWrapRev<FuncWrap, __m256d>;

#if AADC_512
    typedef __m512d mm512Type;
#else
  // fall back to avx256 since these functions will never be called
    typedef __m256d mm512Type;
#endif
    auto call_512_fwd_ptr = &callFuncWrapFwd<FuncWrap, mm512Type>;
    auto call_512_rev_ptr = &callFuncWrapRev<FuncWrap, mm512Type>;

    CAAD_ExtFunction(
      (void*)obj_ptr, (void*)call_256_fwd_ptr, (void*)call_256_rev_ptr
      , (void*)call_512_fwd_ptr, (void*)call_512_rev_ptr
    );
}

};


#endif // AADC_EXT_H
