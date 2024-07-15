#pragma once

#include "aadc/idouble.h"
#include "aadc/ibool.h"
#include "aadc/iint.h"

namespace aadc {



template<class T>
class TypeTraits {
public:
	typedef bool mbool;
	typedef int  mint;
	typedef double mdouble;
};

template<>
class TypeTraits<idouble> {
public:
	typedef ibool mbool;
	typedef iint  mint;
	typedef idouble mdouble;
};

}
