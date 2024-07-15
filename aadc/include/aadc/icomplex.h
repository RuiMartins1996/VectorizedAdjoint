#ifndef CAAD_COMPLEXIDOUBLE_H
#define CAAD_COMPLEXIDOUBLE_H

#include <aadc/idouble.h>
#include <aadc/ibool.h>
#include <complex>

namespace std {

template<>
class complex<idouble>  {
    typedef idouble T;
public:
    complex() {re=0.; im=0.;}
    complex(const double val) : re(val) {im=0;}
    complex(const double r, const double i) : re(r), im(i) {}
    complex(const idouble val) : re(val) {im=0;}
    complex(const idouble r, const idouble i) : re(r), im(i) {}

   complex<T>& operator = (const complex<T>& other) {
        re=other.real();
        im=other.imag();
        return *this;
    }

    complex<T>& operator= (const complex<double>& other) {
        re = other.real();
        im = other.imag();
        return *this;
    }
    
    complex& operator=(double x) {
        re=idouble(x);  //??
        im=0;
        return *this; 
    }

    
    
    const T& real() const { return re; }
    const T& imag() const { return im; }
  
    complex& operator+=(const T& other) { re+=other; return *this; }
    complex& operator-=(const T& other) { re-=other; return *this; }
    complex& operator*=(const T& other) { re*=other; im*=other; return *this; }
    complex& operator/=(const T& other) { re/=other; im/=other; return *this; }

    template<class X>
    complex& operator+=(const std::complex<X>& other) { re+=other.real(); im+=other.imag(); return *this; }
    template<class X>
    complex& operator-=(const std::complex<X>& other) { re-=other.real(); im-=other.imag(); return *this; }
    template<class X>
    complex& operator*=(const std::complex<X>& other) { 
        X re_=re*other.real() - im*other.imag(); 
        im=re*other.imag() + im*other.real();
        re=re_;
        return *this; 
    }
    template<class X>
    complex& operator/=(const std::complex<X>& other) { 
        idouble norm=other.real()*other.real() + other.imag()*other.imag();
        X re_=(re*other.real() + im*other.imag())/norm; 
        im=(-re*other.imag() + im*other.real())/norm;
        re=re_;
        return *this; 
    }

private:
    T re, im;
};

inline std::complex<idouble> iIf(const ibool& cond, const std::complex<idouble>& a, const std::complex<idouble>& b) {
    return std::complex<idouble>(
        iIf(cond, a.real(), b.real())
        , iIf(cond, a.imag(), b.imag())
    ); 
}

inline std::complex<idouble> operator+(const std::complex<idouble>& val ) { 
    return val; 
}
inline std::complex<idouble> operator-(const std::complex<idouble>& val ) { 
    return std::complex<idouble>(-val.real(), -val.imag());
}

inline std::complex<idouble> operator+(const std::complex<idouble>& lhs, const std::complex<idouble>& rhs) {
      return complex<idouble>(lhs.real()+rhs.real(),lhs.imag()+rhs.imag()); 
}
inline std::complex<idouble> operator+(const std::complex<idouble>& lhs, const idouble& rhs) { 
      return complex<idouble>(lhs.real()+rhs,lhs.imag()); 
}
inline std::complex<idouble> operator+(const idouble& lhs, const std::complex<idouble>& rhs) {
      return complex<idouble>(lhs+rhs.real(),rhs.imag()); 
}
inline std::complex<idouble> operator+(const std::complex<idouble>& lhs, const double rhs) {
      return complex<idouble>(lhs.real()+rhs,lhs.imag()); 
}
inline std::complex<idouble> operator+(const double lhs, const std::complex<idouble>& rhs) { 
      return complex<idouble>(lhs+rhs.real(),rhs.imag()); 
}

//--

inline std::complex<idouble> operator-(const std::complex<idouble>& lhs, const std::complex<idouble>& rhs) {
      return complex<idouble>(lhs.real()-rhs.real(),lhs.imag()-rhs.imag()); 
}
inline std::complex<idouble> operator-(const std::complex<idouble>& lhs, const idouble& rhs) { 
      return complex<idouble>(lhs.real()-rhs,lhs.imag());
}
inline std::complex<idouble> operator-(const idouble& lhs, const std::complex<idouble>& rhs) {
    return complex<idouble>(lhs-rhs.real(),-rhs.imag()); 
}
inline std::complex<idouble> operator-(const std::complex<idouble>& lhs, const double rhs) {
      return complex<idouble>(lhs.real()-rhs,lhs.imag()); 
}
inline std::complex<idouble> operator-(const double lhs, const std::complex<idouble>& rhs) { 
      return complex<idouble>(lhs-rhs.real(),-rhs.imag()); 
}
//**
inline std::complex<idouble> operator*(const std::complex<idouble>& lhs, const std::complex<idouble>& rhs) {
        return complex<idouble>(
            lhs.real()*rhs.real() - lhs.imag()*rhs.imag(),  
            lhs.real()*rhs.imag() + lhs.imag()*rhs.real()
        );
}

inline std::complex<idouble> operator*(const std::complex<double>& lhs, const std::complex<idouble>& rhs) {
        return complex<idouble>(
            lhs.real()*rhs.real() - lhs.imag()*rhs.imag(),  
            lhs.real()*rhs.imag() + lhs.imag()*rhs.real()
        );
}

inline std::complex<idouble> operator*(const std::complex<idouble>& lhs, const std::complex<double>& rhs) {
        return complex<idouble>(
            lhs.real()*rhs.real() - lhs.imag()*rhs.imag(),  
            lhs.real()*rhs.imag() + lhs.imag()*rhs.real()
        );
}

inline std::complex<idouble> operator*(const std::complex<idouble>& lhs, const idouble& rhs) { 
        return complex<idouble>(lhs.real()*rhs, lhs.imag()*rhs);
}
inline std::complex<idouble> operator*(const idouble& lhs, const std::complex<idouble>& rhs) {
        return complex<idouble>(rhs.real()*lhs, rhs.imag()*lhs);
}
inline std::complex<idouble> operator*(const std::complex<idouble>& lhs, const double rhs) {
        return complex<idouble>(lhs.real()*rhs, lhs.imag()*rhs);
}
inline std::complex<idouble> operator*(const double lhs, const std::complex<idouble>& rhs) { 
        return complex<idouble>(rhs.real()*lhs, rhs.imag()*lhs);
}
//
inline std::complex<idouble> operator/(const std::complex<idouble>& lhs, const std::complex<idouble>& rhs) {
        idouble norm=1/(rhs.real()*rhs.real() + rhs.imag()*rhs.imag());
        return complex<idouble>(
            (lhs.real()*rhs.real() + lhs.imag()*rhs.imag())*norm,  
            (-lhs.real()*rhs.imag() + lhs.imag()*rhs.real())*norm
        );
}

inline std::complex<idouble> operator/(const std::complex<idouble>& lhs, const std::complex<double>& rhs) {
    idouble norm = 1/(rhs.real()*rhs.real() + rhs.imag()*rhs.imag());
    return complex<idouble>(
        (lhs.real()*rhs.real() + lhs.imag()+rhs.imag())*norm, (-lhs.real()*rhs.imag() + lhs.imag()*rhs.real())*norm
    );
}

inline std::complex<idouble> operator/(const std::complex<idouble>& lhs, const idouble& rhs) { 
        return complex<idouble>(lhs.real()/rhs, lhs.imag()/rhs);
}
inline std::complex<idouble> operator/(const idouble& lhs, const std::complex<idouble>& rhs) {
        idouble norm=1.0/(rhs.real()*rhs.real() + rhs.imag()*rhs.imag());
        return complex<idouble>(lhs*rhs.real()*norm, -lhs*rhs.imag()*norm);
}
inline std::complex<idouble> operator/(const std::complex<idouble>& lhs, const double rhs) {
        return complex<idouble>(lhs.real()/rhs, lhs.imag()/rhs);
}
inline std::complex<idouble> operator/(const double lhs, const std::complex<idouble>& rhs) { 
        idouble norm=1.0/(rhs.real()*rhs.real() + rhs.imag()*rhs.imag());
        return complex<idouble>(lhs*rhs.real()*norm, -lhs*rhs.imag()*norm);
}
//

inline ibool operator==( const complex<idouble>& lhs, const complex<idouble>& rhs) {
    using namespace aadcBoolOps;
    return lhs.real() == rhs.real() && lhs.imag() == rhs.imag(); 
}
inline ibool operator==( const complex<idouble>& lhs, const idouble& rhs) {
    using namespace aadcBoolOps;
    return lhs.real() == rhs && lhs.imag() == 0.; 
}
inline ibool operator==( const idouble& lhs, const complex<idouble>& rhs) { 
    using namespace aadcBoolOps;
    return lhs == rhs.real() && rhs.imag() == 0.; 
}
inline ibool operator!=( const complex<idouble>& lhs, const complex<idouble>& rhs) { 
    using namespace aadcBoolOps;
    return lhs.real() != rhs.real() || lhs.imag() != rhs.imag(); 
}
inline ibool operator!=( const complex<idouble>& lhs, const idouble& rhs) { 
    using namespace aadcBoolOps;
    return lhs.real() != rhs || lhs.imag() != 0.;
}
inline ibool operator!=( const idouble& lhs, const complex<idouble>& rhs) {
    using namespace aadcBoolOps;
    return rhs.imag() != 0. || lhs != rhs.real(); 
}


inline const idouble& real( const std::complex<idouble>& z ) { return z.real(); }
inline const idouble& imag( const std::complex<idouble>& z ) { return z.imag(); }

inline idouble abs( const std::complex<idouble>& z ) { return sqrt(z.real()*z.real() + z.imag()*z.imag());}
inline idouble arg( const std::complex<idouble>& z ) { return atan2(z.imag(), z.real());}

inline idouble norm( const std::complex<idouble>& z ) { return z.real()*z.real() + z.imag()*z.imag(); }
inline std::complex<idouble> conj( const std::complex<idouble>& z ) { return std::complex<idouble>(z.real(), -z.imag()); }

//inline std::complex<idouble> proj( const std::complex<idouble>& z ) {
//    const idouble den = z.real()*z.real() + z.imag()*z.imag() + 1.0;
//    return std::complex<idouble>((2.0 * z.real()) / den, (2.0 * z.imag()) / den);
//}

inline std::complex<idouble> polar(const idouble& r, const idouble& theta = 0.) { 
    //__glibcxx_assert( __rho >= 0 );
    return std::complex<idouble>(r*cos(theta), r*sin(theta));
}

inline std::complex<idouble> exp( const std::complex<idouble>& z ) {return std::polar(exp(z.real()), z.imag());}
inline std::complex<idouble> log( const std::complex<idouble>& z ) {return std::complex<idouble>(log(std::abs(z)), std::arg(z));}
inline std::complex<idouble> log10( const std::complex<idouble>& z ) { return std::log(z) / log(10.0);}
inline std::complex<idouble> sqrt( const std::complex<idouble>& z ) { 
    const idouble& x = z.real();
    const idouble& y = z.imag();
    idouble t = sqrt(2.0 * (std::abs(z) + abs(x)));
    idouble u = t * 0.5;
//    idouble t0 = sqrt(abs(y) * 0.5);
    return iIf(
        (x == 0.), complex<idouble>(u, iIf(y < 0.,  -u, u)), 
        iIf(
            x > 0., complex<idouble>(u, y / t)
            , complex<idouble>(abs(y) / t, iIf(y < 0., -u, u)) 
        )
    );
} 

inline std::complex<idouble> pow( const std::complex<idouble>& x, const std::complex<idouble>& y) { 
    using namespace aadcBoolOps;
    return iIf((x.real() ==0.) && (x.imag() ==0. ) , std::complex<idouble>(), std::exp(y * std::log(x))); 
}

inline std::complex<idouble> pow( const std::complex<idouble>& x, const idouble& y) { 
    using namespace aadcBoolOps;
    complex<idouble> t = std::log(x);
    return iIf ( 
        (x.imag() == 0. && x.real() > 0.) ,  pow(x.real(), y), 
        std::polar(exp(y * t.real()), y * t.imag())
    );     
}
inline std::complex<idouble> pow( const idouble& x, const complex<idouble>& y) { 
    using namespace aadcBoolOps;
    return iIf(x == 0 , 0., std::exp(y * std::log(x)));}

inline std::complex<idouble> pow( const complex<idouble>& x, double y) { return pow(x,idouble(y)); }
inline std::complex<idouble> pow( double x, const complex<idouble>& y) { return pow(idouble(x),y);  }

inline complex<idouble> sin( const complex<idouble>& z ) { 
    const idouble x = z.real();
    const idouble y = z.imag();
    return complex<idouble>(sin(x) * cosh(y), cos(x) * sinh(y));
}
inline complex<idouble> cos( const complex<idouble>& z ) {
    const idouble x = z.real();
    const idouble y = z.imag();
    return  complex<idouble>(cos(x) * cosh(y), -sin(x) * sinh(y));
}
inline complex<idouble> tan( const complex<idouble>& z ) { return std::sin(z) / std::cos(z); }

inline complex<idouble> asinh( const complex<idouble>& z ) { 
    std::complex<idouble> t((z.real() - z.imag()) * (z.real() + z.imag()) + 1.0, 2.0 * z.real() * z.imag());
    t = std::sqrt(t);
    return std::log(t + z);
}

inline complex<idouble> asin( const complex<idouble>& z ) {  
    complex<idouble> t(-z.imag(), z.real());
    t = std::asinh(t);
    return std::complex<idouble>(t.imag(), -t.real());
}
inline complex<idouble> acos( const complex<idouble>& z ) { 
    const std::complex<idouble> t = std::asin(z);
    const idouble pi_2 = 1.5707963267948966192313216916397514L;
    return std::complex<idouble>(pi_2 - t.real(), -t.imag());

}
inline complex<idouble> atan( const complex<idouble>& z ) { 
    const idouble r2 = z.real() * z.real();
    const idouble x = 1.0 - r2 - z.imag() * z.imag();
    idouble num = z.imag() + 1.0;
    idouble den = z.imag() - 1.0;
    num = r2 + num * num;
    den = r2 + den * den;
    return std::complex<idouble>(0.5 * atan2(2.0 * z.real(), x), 0.25 * log(num / den));
}

inline complex<idouble> sinh( const complex<idouble>& z ) { return 0.5 * (exp(z)-exp(-z)); }
inline complex<idouble> cosh( const complex<idouble>& z ) { return 0.5 * (exp(z)+exp(-z)); }
inline complex<idouble> tanh( const complex<idouble>& z ) { return std::sinh(z) / std::cosh(z); }

inline complex<idouble> acosh( const complex<idouble>& z ) { 
    return 2.0 * std::log(std::sqrt(0.5 * (z + 1.0)) + std::sqrt(0.5 * (z - 1.0)));
}
inline complex<idouble> atanh( const complex<idouble>& z ) { 
    const idouble i2 = z.imag() * z.imag();
    const idouble x = 1.0 - i2 - z.real() * z.real();
    idouble num = 1.0 + z.real();
    idouble den = 1.0 - z.real();
    num = i2 + num * num;
    den = i2 + den * den;
    return std::complex<idouble>(0.25 * (log(num) - log(den)) , 0.5 * atan2(2.0 * z.imag(), x));
}


};
#endif // CAAD_COMPLEXIDOUBLE_H