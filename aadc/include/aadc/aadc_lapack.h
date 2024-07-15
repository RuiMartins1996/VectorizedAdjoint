#pragma once
#include <aadc/idouble.h>
#include <aadc/iint.h>
#include <aadc/aadc_ext.h>
#include "lapack.h"

//dposv
//https://www.nag.com/numeric/fl/nagdoc_latest/html/f07/f07faf.html#uplo


///////////////////////////////////////////////////////////////////////
//                   DGETRS
//  https://www.nag.com/numeric/fl/nagdoc_latest/html/f07/f07aef.html
///////////////////////////////////////////////////////////////////////
class AadcDgetrsWrapper : public aadc::ConstStateExtFunc {
public:
    AadcDgetrsWrapper(//bool first_avx_only, 
        char *trans, int *n,
        int *nrhs, const idouble *A, int *ldA, 
        const iint *ipiv, idouble *B, int *ldB, int *info, idouble *B_out
    ) : trans_(*trans), n_(*n), nrhs_(*nrhs), ldA_(*ldA), ldB_(*ldB), info_(*info) 
    {
            A_indx.reserve(n_ * n_);
            B_out_indx.reserve(n_ * nrhs_);
            ipiv_indx.reserve(n_);

            bool diff_check = false;
            for (int i=0; i < n_ * n_; i++) {
                diff_check = diff_check || CAAD_iVarIsDiff(&(A[i].val));
                A_indx.emplace_back(ExtVarIndex(A[i], true));
            }
            for (int i=0; i < n_ * nrhs_; i++) {
                diff_check = diff_check || CAAD_iVarIsDiff(&(B[i].val));
                B_indx.emplace_back(B[i], false);
            }
            for (int i=0; i < n_ * nrhs_; i++) {
                B_out_indx.emplace_back(ExtVarIndex(B_out[i], true, diff_check));
            }
            for (int i=0; i < n_; i++) ipiv_indx.emplace_back(ExtIntVarIndex(ipiv[i], true));
        std::vector<int> to_int_ipiv(n_);
        for (int i=0; i<n_; i++) to_int_ipiv[i]=int(ipiv[i].val);
        dgetrs_(&trans_, &n_, &nrhs_, (double *)A, &ldA_, &to_int_ipiv.front(), (double *)B_out, &ldB_, &info_);
        for(int i=0; i < n_ * nrhs_; i++) B_out[i].val = B[i].val;
    }
    
    template<typename mmType>
    void forward(mmType* v) const {
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(B_indx.size());
        std::vector<double> argA(A_indx.size());
        std::vector<int> arg_ipiv(ipiv_indx.size()); 
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(v[B_indx[i]])[avx_i];
            for (int i=0; i < argA.size(); i++) argA[i] = aadc::toDblPtr(v[A_indx[i]])[avx_i];
            for (int i=0; i < ipiv_indx.size(); i++) arg_ipiv[i] = ((int64_t*)(aadc::toDblPtr(v[ipiv_indx[i]])))[avx_i];
            int info_aux; // forward const
            dgetrs_(&trans_, &n_, &nrhs_,  &argA.front(), &ldA_, &arg_ipiv.front(), &res.front(), &ldB_, &info_aux);
            for (int i=0; i<res.size(); i++)  aadc::toDblPtr(v[B_out_indx[i]])[avx_i] = res[i];
        }
    }

    template<class mmType>
    void reverse(const mmType *v, mmType *d) const {
        char trans_inv = trans_ == 'N' ? 'T' : 'N';
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(B_indx.size());
        std::vector<double> argA(A_indx.size());
        std::vector<int> arg_ipiv(ipiv_indx.size()); 
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            std::vector<double> dA(A_indx.size(), 0.);
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(d[B_out_indx[i]])[avx_i];
            for (int i=0; i < argA.size(); i++) argA[i] = aadc::toDblPtr(v[A_indx[i]])[avx_i];
            for (int i=0; i < ipiv_indx.size(); i++) arg_ipiv[i] = ((int64_t*)(aadc::toDblPtr(v[ipiv_indx[i]])))[avx_i];
            int info_aux; 

            //C=A^{-1}B. here we calculate dB. If trans != N, then C=A^{-T} B
            dgetrs_(&trans_inv, &n_, &nrhs_, &argA.front(), &ldA_, &arg_ipiv.front(), &res.front(), &ldB_, &info_aux);
            for (int i=0; i<res.size(); i++) aadc::toDblPtr(d[B_indx[i]])[avx_i] += res[i];

            // A = PLU. Here we calculate A_bar
            
            int i=0, j=0;
            for (int m=0; m<argA.size(); m++) {
                for (int k=0; k<nrhs_; k++) {
                    dA[m] += -res[i + k*n_] * aadc::toDblPtr(v[B_out_indx[j + k*n_]])[avx_i];
                }
                i++;
                if (i==n_) {
                    j++;
                    i=0;
                }
            }

            // if trans != N, then we transpose dA
            if (trans_inv == 'N') {
                for (int i=0; i< n_; i++) {
                    for (int j=0; j<i; j++) {
                        double tmp = dA[i*n_+ j];
                        dA[i*n_+j]= dA[j*n_ + i];
                        dA[j*n_+i] = tmp;
                    }
                }
            }

            // compute matrix P through the ipiv, i.e. P_{i,j}=1 iff j=final_indx[i] 
            std::vector<int> final_indx(ipiv_indx.size()); 
            for (int i=0; i<ipiv_indx.size(); i++) final_indx[i]=i;
            for (int i=0; i<ipiv_indx.size(); i++) {
                int tmp = final_indx[i];
                int k = arg_ipiv[i]-1;
                final_indx[i] = final_indx[k];
                final_indx[k] = tmp;
            }
            
            // L_bar =P^t A_bar U^t
            for (int i=0; i< n_; i++) {
                for (int j=0; j<i; j++) {
                    int dpos = j*n_ + i;
                    for (int k=j; k<n_; k++) {
                        int int_pos = final_indx[i] + n_*k;
                        int int_pos2 = k*n_ + j;
                        aadc::toDblPtr(d[A_indx[dpos]])[avx_i] += dA[int_pos] * argA[int_pos2];
                    }
                }
            }

            // U_bar = L^t P^t A_bar
            for (int i=0; i< n_; i++) {
                for (int j=i; j<n_; j++) {
                    int dpos = i + n_*j;
                    for (int k=i; k< n_; k++) {
                        int int_pos = i*n_ + k;
                        int int_pos2 = j*n_ + final_indx[k];
                        aadc::toDblPtr(d[A_indx[dpos]])[avx_i] += (k==i ? 1 : argA[int_pos]) * dA[int_pos2];
                    }
                }
            }
        }
    }

private:
    std::vector<ExtVarIndex> A_indx, B_indx, B_out_indx;
    std::vector<ExtIntVarIndex> ipiv_indx;
    char trans_;
    int n_;
    int nrhs_; 
    int ldA_;
    int ldB_;
    int trans_len_;
    int info_;
};

inline void dgetrs_(//bool first_avx_only,
    char *trans, int *n,
    int *nrhs, const idouble *A, int *ldA,
    const iint *ipiv, idouble *B, int *ldB, int *info  
) {
    if (!idouble::recording) {
        std::vector<int> to_int_ipiv(*n);
        for (int i=0; i<*n; i++) to_int_ipiv[i]=int(ipiv[i].val);
        dgetrs_(trans, n, nrhs, (const double *)A, ldA, (int *) (&to_int_ipiv.front()), (double *)B, ldB, info);
        return;
    }
    std::vector<idouble> B_out((*n) * (*nrhs));
    auto wrapper(
        std::make_shared<AadcDgetrsWrapper>(
            trans, n, nrhs, A, ldA, ipiv, B, ldB, info, &B_out.front() 
        )
    );
    aadc::addConstStateExtFunction(wrapper);
    for (int i=0; i < B_out.size(); i++) B[i] = B_out[i];
}
 
////////////////////////////////////////////////////////////////////////////////
//             DGETRF
// https://www.nag.com/numeric/fl/nagdoc_latest/html/f07/f07adf.html
///////////////////////////////////////////////////////////////////////////////
class AadcDgetrfWrapper : public aadc::ConstStateExtFunc {
public:
    AadcDgetrfWrapper(
        int *m, int *n,  idouble *A, int *ldA, iint *ipiv, int *info, idouble *A_out
        ) : m_(*m), n_(*n), ldA_(*ldA), info_(*info) 
    {
        int ipiv_size = std::min(m_, n_);
        std::vector<int> to_int_ipiv(ipiv_size);
        for (int i=0; i < ipiv_size; i++) to_int_ipiv[i] = int(ipiv[i].val);
        
        bool out_diff = false;
        A_indx.reserve(m_ * n_);
        A_diff.reserve(m_ * n_);
        A_out_indx.reserve(m_ * n_);
        ipiv_indx.reserve(ipiv_size);
        
        for (int i=0; i<m_ * n_; i++) {
            A_diff.emplace_back(CAAD_iVarIsDiff(&(A[i].val)));
            out_diff = out_diff || A_diff.back();
            A_indx.emplace_back(ExtVarIndex(A[i], true));
        }
        for (int i=0; i<m_ * n_; i++) A_out_indx.emplace_back(ExtVarIndex(A_out[i], true, out_diff));
        for (int i=0; i < ipiv_size; i++) ipiv_indx.emplace_back(ExtIntVarIndex(ipiv[i], true));
    
        
        dgetrf_(&m_, &n_, (double*)A, &ldA_, &to_int_ipiv.front(), info);
        for (int i=0; i<to_int_ipiv.size(); i++) ipiv[i].val = to_int_ipiv[i];
        for (int i=0; i<m_ * n_; i++) A_out[i].val = A[i].val;
    } 

    template<typename mmType>
    void forward(mmType* v) const  {
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(A_indx.size());
        std::vector<int> res_ipiv(ipiv_indx.size()); 
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(v[A_indx[i]])[avx_i];
            int info;
            dgetrf_(&m_, &n_,  &*res.begin(), &ldA_,  &*res_ipiv.begin(), &info);
            for (int i=0; i < res.size(); i++) aadc::toDblPtr(v[A_out_indx[i]])[avx_i] = res[i];
            for (int i=0; i < res_ipiv.size(); i++) {
                ((int64_t*)(aadc::toDblPtr(v[ipiv_indx[i]])))[avx_i] = res_ipiv[i];
            }
        }
    }

    template<class mmType>
    void reverse(const mmType *v, mmType *d) const {
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(A_indx.size());
        double h = 1e-6;
        std::vector<int> res_ipiv(ipiv_indx.size()), res_ipiv_bump(ipiv_indx.size()); 
        int info;
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(v[A_indx[i]])[avx_i];
            std::vector<double> base(res);
            dgetrf_(&m_, &n_,  &*res.begin(), &ldA_,  &*res_ipiv.begin(), &info);
            for (int j=0; j<A_indx.size(); j++) {
                if (!A_diff[j]) continue;
                std::vector<double> res_bump(base);
                res_bump[j] += h;
                dgetrf_(&m_, &n_, &*res_bump.begin(), &ldA_, &*res_ipiv_bump.begin(), &info);

                //bool check = std::equal(res_ipiv.begin(), res_ipiv.end(), res_ipiv_bump.begin());
                bool check = false;
                for (int i=0; i<res_ipiv_bump.size(); i++) if (res_ipiv_bump[i] != res_ipiv[i]) check = true;
                if (check)  {
                    for (int i=0; i<res_ipiv_bump.size(); i++) std::cout << res_ipiv[i] << " " << res_ipiv_bump[i] << std::endl;
                    throw std::runtime_error("Ipiv changed after bump during the dgetrf reverse");
                }

                /*
                std::vector<double> res_minus_bump(base);
                res_minus_bump[j] -= h;
                dgetrf_(&m_, &n_, &*res_minus_bump.begin(), &ldA_, &*res_ipiv_bump.begin(), &info);

                //bool check = std::equal(res_ipiv.begin(), res_ipiv.end(), res_ipiv_bump.begin());
                check = false;
                for (int i=0; i<res_ipiv_bump.size(); i++) if (res_ipiv_bump[i] != res_ipiv[i]) check = true;
                if (check)  {
                    for (int i=0; i<res_ipiv_bump.size(); i++) std::cout << res_ipiv[i] << " " << res_ipiv_bump[i] << std::endl;
                    throw std::runtime_error("Ipiv changed after bump during the dgetrf reverse");
                }

                // +2h
                std::vector<double> res_upup_bump(base);
                res_upup_bump[j] += 2*h;
                dgetrf_(&m_, &n_, &*res_upup_bump.begin(), &ldA_, &*res_ipiv_bump.begin(), &info);

                //bool check = std::equal(res_ipiv.begin(), res_ipiv.end(), res_ipiv_bump.begin());
                check = false;
                for (int i=0; i<res_ipiv_bump.size(); i++) if (res_ipiv_bump[i] != res_ipiv[i]) check = true;
                if (check)  {
                    for (int i=0; i<res_ipiv_bump.size(); i++) std::cout << res_ipiv[i] << " " << res_ipiv_bump[i] << std::endl;
                    throw std::runtime_error("Ipiv changed after bump during the dgetrf reverse");
                }

                //
                std::vector<double> res_minusminus_bump(base);
                res_minusminus_bump[j] -= 2*h;
                dgetrf_(&m_, &n_, &*res_minusminus_bump.begin(), &ldA_, &*res_ipiv_bump.begin(), &info);

                //bool check = std::equal(res_ipiv.begin(), res_ipiv.end(), res_ipiv_bump.begin());
                check = false;
                for (int i=0; i<res_ipiv_bump.size(); i++) if (res_ipiv_bump[i] != res_ipiv[i]) check = true;
                if (check)  {
                    for (int i=0; i<res_ipiv_bump.size(); i++) std::cout << res_ipiv[i] << " " << res_ipiv_bump[i] << std::endl;
                    throw std::runtime_error("Ipiv changed after bump during the dgetrf reverse");
                }
                */
                for (int l=0; l<A_indx.size(); l++) {
                    aadc::toDblPtr(d[A_indx[j]])[avx_i] += (res_bump[l] - res[l]) / 1 / h * aadc::toDblPtr(d[A_out_indx[l]])[avx_i];
                        //(-res_upup_bump[l] + 8*res_bump[l] - 8 * res_minus_bump[l] + res_minusminus_bump[l]) / 12 / h 
                        //* aadc::toDblPtr(d[A_out_indx[l]])[avx_i]
                    //; 
                }
            }
        }
    }

private:
    std::vector<ExtVarIndex> A_indx, A_out_indx;
    std::vector<ExtIntVarIndex> ipiv_indx;
    std::vector<bool> A_diff;
    int m_,n_,ldA_;
    int info_;
};

inline void dgetrf_(int *m, int *n,  idouble *A, int *ldA, iint *ipiv, int *info){
    // TODO: ADD case if recording isnt't taking place!!!
    if (!idouble::recording) {
        std::vector<int> to_int_ipiv(*n);
        for (int i=0; i<*n; i++) to_int_ipiv[i]=int(ipiv[i].val);
        dgetrf_(m, n, (double*)A, ldA,  (int *) (&to_int_ipiv.front()), info);
        for (int i=0; i<*n; i++) ipiv[i]=to_int_ipiv[i];
        return;
    }
    
    std::vector<idouble> A_out((*m) * (*n));
    auto wrapper(
        std::make_shared<AadcDgetrfWrapper>(
         m, n, A, ldA, ipiv, info, &A_out.front()
    ));
    aadc::addConstStateExtFunction(wrapper);
    for (int i=0; i < A_out.size(); i++) A[i] = A_out[i];
}


////////////////////////////////////////////////////////////////////////////////
//             linSolver
// Similar to Lapack dgesv routine, but A is fixed in linSolver
///////////////////////////////////////////////////////////////////////////////

inline void linSolver(
    const char *trans, const int *n,
    const int *nrhs,  double *A, const int *ldA, 
    double *B, const int *ldB, int *info
){    
    std::vector<double> A_tmp((*n) * (*n));
    for (int i=0; i < A_tmp.size(); i++) A_tmp[i] = A[i];
    
    std::vector<int> ipiv(*n);
    dgetrf_(n, n, &A_tmp.front(), ldA, &ipiv.front(), info);
    dgetrs_(trans, n, nrhs, &A_tmp.front(), ldA, &ipiv.front(), B, ldB, info);
    //dgetrf_(n, n, A, ldA, ipiv, info);
    //dgetrs_(trans, n, nrhs, A, ldA, ipiv, B, ldB, info);
}

class AadcLinSolverWrapper : public aadc::ConstStateExtFunc {
public:
    AadcLinSolverWrapper(//bool first_avx_only, 
        const char *trans, int *n,
        int *nrhs, const idouble *A, const int *ldA, 
        const iint *ipiv, idouble *B, const int *ldB, int *info, idouble *B_out,  idouble *A_out
    ) : trans_(*trans), n_(*n), nrhs_(*nrhs), ldA_(*ldA), ldB_(*ldB), info_(*info) 
    {
        A_indx.reserve(n_ * n_);
        A_out_indx.reserve(n_ * n_);
        B_out_indx.reserve(n_ * nrhs_);
        B_indx.reserve(n_ * nrhs_);
        ipiv_indx.reserve(n_);

        for (int i=0; i < n_ * n_; i++) A_indx.emplace_back(ExtVarIndex(A[i], false));
        for (int i=0; i < n_ * nrhs_; i++) B_indx.emplace_back(ExtVarIndex(B[i], false));
        for (int i=0; i < n_ * nrhs_; i++) B_out_indx.emplace_back(ExtVarIndex(B_out[i], true));
        for (int i=0; i < n_; i++) ipiv_indx.emplace_back(ExtIntVarIndex(ipiv[i], true));
        for (int i=0; i<n_ * n_; i++) A_out_indx.emplace_back(ExtVarIndex(A_out[i], true));

        linSolver(&trans_, &n_, &nrhs_, (double*)A, &ldA_, (double *)B_out, &ldB_, &info_);
        for (int i=0; i<n_ * n_; i++) A_out[i].val = A[i].val;
        for (int i=0; i < n_ * nrhs_; i++) B_out[i].val = B[i].val;
    }
    
    template<typename mmType>
    void forward(mmType* v) const  {
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(B_indx.size());
        std::vector<double> argA(A_indx.size());
        std::vector<int> arg_ipiv(ipiv_indx.size()); 
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(v[B_indx[i]])[avx_i];
            for (int i=0; i < argA.size(); i++) argA[i] = aadc::toDblPtr(v[A_indx[i]])[avx_i];
            int info_aux; // forward const

            dgetrf_(&n_, &n_,  &*argA.begin(), &ldA_,  &*arg_ipiv.begin(), &info_aux);
            for (int i=0; i < argA.size(); i++) aadc::toDblPtr(v[A_out_indx[i]])[avx_i] = argA[i];
            for (int i=0; i < arg_ipiv.size(); i++) {
                ((int64_t*)(aadc::toDblPtr(v[ipiv_indx[i]])))[avx_i] = arg_ipiv[i];
            }

            dgetrs_(&trans_, &n_, &nrhs_,  &argA.front(), &ldA_, &arg_ipiv.front(), &res.front(), &ldB_, &info_aux);
            for (int i=0; i<res.size(); i++)  aadc::toDblPtr(v[B_out_indx[i]])[avx_i] = res[i];
        }
    }

    template<class mmType>
    void reverse(const mmType *v, mmType *d) const {
        char trans_inv = trans_ == 'N' ? 'T' : 'N';
        const int AVX_length = aadc::mmSize<mmType>();
        std::vector<double> res(B_indx.size());
        std::vector<double> argA(A_indx.size());
        std::vector<int> arg_ipiv(ipiv_indx.size()); 
        for (int avx_i=0; avx_i<AVX_length; avx_i++) {
            for (int i=0; i < res.size(); i++) res[i] = aadc::toDblPtr(d[B_out_indx[i]])[avx_i];
            for (int i=0; i < argA.size(); i++) argA[i] = aadc::toDblPtr(v[A_out_indx[i]])[avx_i];
            for (int i=0; i < ipiv_indx.size(); i++) arg_ipiv[i] = ((int64_t*)(aadc::toDblPtr(v[ipiv_indx[i]])))[avx_i];
            int info_aux; // forward const

            dgetrs_(&trans_inv, &n_, &nrhs_, &argA.front(), &ldA_, &arg_ipiv.front(), &res.front(), &ldB_, &info_aux);
            for (int i=0; i<res.size(); i++) aadc::toDblPtr(d[B_indx[i]])[avx_i] += res[i];

            if (trans_inv != 'N') {
                int i=0,j=0;
                for (int m=0; m<argA.size(); m++) {
                    for (int k=0; k<nrhs_; k++) {
                        aadc::toDblPtr(d[A_indx[m]])[avx_i] += -res[i + k*n_] * aadc::toDblPtr(v[B_out_indx[j+ k*n_]])[avx_i];
                    }
                    i++;
                    if (i==n_) {
                        j++;
                        i=0;
                    }
                }
            } else {
                int i=0,j=0;
                for (int m=0; m<argA.size(); m++) {
                    for (int k=0; k<nrhs_; k++) {
                        aadc::toDblPtr(d[A_indx[m]])[avx_i] += -res[i + k*n_] * aadc::toDblPtr(v[B_out_indx[j + k*n_]])[avx_i];
                    }
                    j++;
                    if (j==n_) {
                        i++;
                        j=0;
                    }
                }
            }                    
        }
    }

private:
    std::vector<ExtVarIndex> A_indx, B_indx, B_out_indx, A_out_indx;
    std::vector<ExtIntVarIndex> ipiv_indx;
    char trans_;
    int n_;
    int nrhs_; 
    int ldA_;
    int ldB_;
    int trans_len_;
    int info_;
};

inline void linSolver(//bool first_avx_only,
    char *trans, int *n,
    int *nrhs, idouble *A, int *ldA,
    idouble *B, int *ldB, int *info  
) {
    if (!idouble::recording) {
        linSolver(trans, n, nrhs, (double*)A, ldA, (double *)B, ldB, info);
        return;
    }    
  
    std::vector<idouble> A_tmp((*n) * (*n));
    for (int i=0; i < A_tmp.size(); i++) A_tmp[i] = A[i];

    std::vector<idouble> B_out((*n) * (*nrhs));
    std::vector<idouble> A_out((*n) * (*n));
    std::vector<iint> ipiv(*n);

    auto wrapper(
        std::make_shared<AadcLinSolverWrapper>(
            trans, n, nrhs, &A_tmp.front(), ldA, &ipiv.front(), B, ldB, info, &B_out.front(), &A_out.front() 
        )
    );
    aadc::addConstStateExtFunction(wrapper);
    for (int i=0; i < B_out.size(); i++) B[i] = B_out[i];
    //for (int i=0; i < A_out.size(); i++) A_tmp[i] = A_out[i];
    //for (int i=0; i < A_out.size(); i++) A[i]= A_out[i];
}