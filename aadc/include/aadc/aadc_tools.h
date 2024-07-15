#pragma once

#include <immintrin.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "tools/mmvector.h"

namespace aadc {

	template<typename mmType>
	inline uint64_t mmSize() {
		return sizeof(mmType) / sizeof(double);
	}

	inline double* toDblPtr(__m256d& v)
	{
		return reinterpret_cast<double*>(&v);
	}

	inline const double* toDblPtr(const __m256d& v)
	{
		return reinterpret_cast<const double*>(&v);
	}
#if AADC_512
	inline double* toDblPtr(__m512d& v)
	{
		return reinterpret_cast<double*>(&v);
	}
	inline const double* toDblPtr(const __m512d& v)
	{
		return reinterpret_cast<const double*>(&v);
	}
#endif
	inline double mmSum(const __m256d& v) {
		const double* vd(toDblPtr(v));
		return vd[0] + vd[1] + vd[2] + vd[3];
	}

    template<typename mmType> inline mmType  mmZero()             { return mmType(0.0);         }
    template<>                inline __m256d mmZero<__m256d>()    { return _mm256_setzero_pd(); }
    template<typename mmType> inline mmType  mmSetConst(double v) { return _mm256_set1_pd(v);   }
	template<> inline double  mmSetConst<double>(double v) { return v;   }

    inline  double mmAdd (const double a, const double b) {
        return a+b;
    }

	inline  __m256d mmAdd (const __m256d a, const __m256d b) {
        return _mm256_add_pd(a,b);
    }
    inline __m256d mmSub (const __m256d a, const __m256d b) {
        return _mm256_sub_pd(a,b);
    }
	inline  __m256d mmMul (const __m256d a, const __m256d b) {
		return _mm256_mul_pd(a,b);
	}
	inline __m256d  mmDiv (const __m256d a, const __m256d b) {
		return _mm256_div_pd(a,b);
	}
#if AADC_512
	inline double mmSum(const __m512d& v) {
		const double* vd(toDblPtr(v));
		return vd[0] + vd[1] + vd[2] + vd[3]
			+ vd[4] + vd[5] + vd[6] + vd[7]
		;
	}

    template<> inline __m512d mmZero<__m512d>    ()         { return _mm512_setzero_pd(); }
    template<> inline __m512d mmSetConst<__m512d>(double v) { return _mm512_set1_pd(v);   }

    inline __m512d mmSub (const __m512d a, const __m512d b) {
        return _mm512_sub_pd(a,b);
    }
    inline __m512d mmAdd (const __m512d a, const __m512d b) {
        return _mm512_add_pd(a,b);
    }
	inline __m512d mmMul (const __m512d a, const __m512d b) {
		return _mm512_mul_pd(a,b);
	}
	inline __m512d mmDiv (const __m512d a, const __m512d b) {
		return _mm512_div_pd(a,b);
	}
#endif

	template<typename mmType>
	mmType mmMax (const mmType& a,const mmType& b) {
		mmType x;
		auto _a = toDblPtr(a);
		auto _b = toDblPtr(b);
		auto _x = toDblPtr(x);
		for (int i=0; i<aadc::mmSize<mmType>(); i++) {
			_x[i]=std::max(_a[i],_b[i]);
		}
		return x;
	}

	// using template as there is no intrinsic _mm256_abs_pd 
	template<typename mmType>
	mmType mmFabs(const mmType& a) {
		mmType x;
		auto _a = toDblPtr(a);
		auto _x = toDblPtr(x);
		for (int i = 0; i < aadc::mmSize<mmType>(); i++) {
			_x[i] = fabs(_a[i]);
		}
		return x;
	}

	template<typename mmType>
	mmType mmIndicatorFunction(const mmType& a) {
		mmType x;
		auto _a = toDblPtr(a);
		auto _x = toDblPtr(x);
		for (int i=0; i<aadc::mmSize<mmType>(); i++) {
			_x[i]=(_a[i]>0.) ? 1 : 0;
		}
		return x;
	}

};
