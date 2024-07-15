#pragma once

#include <aadc/aadc.h>

extern "C" AADC_API void CAAD_MatrixVectorProduct(
    const double* mat, const double* vec, double* res,
    uint64_t num_r, uint64_t num_c
);
namespace aadc {


template<class mdouble>
class Matrix {
public:
    Matrix(uint64_t nr, uint64_t nc)
        : m_nr(nr), m_nc(nc)
        , m_data(m_nr * m_nc)
    {}
    Matrix(uint64_t nr, uint64_t nc, mdouble val)
        : m_nr(nr), m_nc(nc)
        , m_data(m_nr * m_nc, val)
    {}
    Matrix() : Matrix(0,0) {}
public:
    Matrix<mdouble>& operator = (const Matrix& other) {
        m_nr = other.m_nr;
        m_nc = other.m_nc;
        m_data = other.m_data;
        return *this;
    }

    mdouble& operator ()(const uint64_t r, const uint64_t c) {
        return m_data[r * m_nc + c];
    }
    const mdouble& operator ()(const uint64_t r, const uint64_t c) const {
        return m_data[r * m_nc + c];
    }

    mdouble* operator [](const uint64_t r) {
        return &m_data[r * m_nc];
    }
    const mdouble* operator [](const uint64_t r) const {
        return &m_data[r * m_nc];
    }

    void resize(uint64_t nr, uint64_t nc) {
        m_nr = nr; m_nc = nc;
        m_data.resize(nr*nc);
    }

    uint64_t rows() const { return m_nr; }
    uint64_t cols() const { return m_nc; }

private:
    uint64_t m_nr, m_nc;
    std::vector<mdouble> m_data;
};



// typedef Matrix<idouble>std::vector<std::vector<idouble> > iMatrix;
typedef Matrix<idouble> iMatrix;
//typedef std::vector<std::vector<double> > ScalarMatrix;
typedef Matrix<double> ScalarMatrix;
template <class mmType>
using AVXMatrix = std::vector<mmVector<mmType> >;
//typedef std::vector<std::vector<aadc::AADCArgument> > MatrixArg;
//typedef std::vector<std::vector<aadc::AADCScalarArgument> > ScalarMatrixArg;
//typedef std::vector<std::vector<aadc::AADCResult> > MatrixResult;

typedef Matrix<aadc::AADCArgument> MatrixArg;
typedef Matrix<aadc::AADCScalarArgument> ScalarMatrixArg;
typedef Matrix<aadc::AADCResult> MatrixResult;

typedef std::vector<idouble> iVector;
typedef std::vector<double> ScalarVector;
template <class mmType>
using AVXVector = mmVector<mmType>;
typedef std::vector<aadc::AADCArgument> VectorArg;
typedef std::vector<aadc::AADCResult> VectorRes;

typedef std::vector<aadc::AADCScalarArgument> ScalarVectorArg;


template<class Matrix>
void initNNMatrix( Matrix& m, int size_row, int size_col) {
    m.resize(size_row);
    for (int ri = 0; ri < size_row; ++ri) {
        m[ri].resize(size_col);
        for (int ci = 0;ci < size_col; ++ci) {
            m[ri][ci] = double(std::rand()) / RAND_MAX - 0.5;
        }
    }
}

template<class Vector>
void initNNVector( Vector& m, int size) {
    m.resize(size);
    for (int ri = 0; ri < size; ++ri) {
        m[ri] = double(std::rand()) / RAND_MAX - 0.5;
    }
}

template<class mmType>
void initAVXMatrix( AVXMatrix<mmType>& m, int size_row, int size_col) {
    m.resize(size_row);
    for (int ri = 0; ri < size_row; ++ri) {
        m[ri].resize(size_col);
        for (int ci = 0;ci < size_col; ++ci) {
            m[ri][ci] = mmSetConst<mmType>(0.1);
        }
    }
}


inline void markMatrixAsInput(MatrixArg& mi, iMatrix& m, const bool diff) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, m.rows()*m.cols()));

    mi.resize(m.rows(), m.cols());
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            m[ri][ci].initVarAtAddress(addr_start);++addr_start;
            if (diff) {
                mi[ri][ci] = m[ri][ci].markAsInput();
            } else mi[ri][ci] = m[ri][ci].markAsInputNoDiff();
        }
    }
}


inline void markMatrixAsOutput(MatrixResult& mi, iMatrix& m) {
    mi.resize(m.rows(), m.cols());
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            mi[ri][ci] = m[ri][ci].markAsOutput();
        }
    }
}


inline iMatrix transpose(iMatrix& m) {
    iMatrix matrix; 
    matrix.resize(m.rows(), m.cols());
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            matrix[ri][ci] = m[ci][ri];
        }
    }
    return matrix;
}


inline void markScalarMatrixAsInput(ScalarMatrixArg& mi, iMatrix& m) {
    int avx_size = CAAD_GetAVXSize();
    int rem(m.rows() % avx_size);
    int avx_r = (m.rows() / avx_size) + (rem ? 1 : 0);
    int nr = m.rows(), nc = m.cols();

//    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, m.size()*m[0].size()));
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, avx_size*avx_r*nc));

    mi.resize(m.rows(), m.cols());
    for (int avx_ri = 0; avx_ri < avx_r; ++avx_ri) {
        for (int avx_i = 0; avx_i < avx_size; ++avx_i) {
            int ri = avx_ri * avx_size + avx_i;
            if (ri < nr) {
                for (int ci = 0;ci < nc; ++ci) {
                    uint64_t addr(addr_start + avx_ri*avx_size*nc+avx_i+ci*avx_size);
//                    std::cout << "mark at " << ri << " " << ci << " " << (uint64_t)(&m[ri][ci]) << " " << addr << std::endl;
                    m[ri][ci].initScalarVarAtAddress(addr_start + avx_ri*avx_size*nc+avx_i+ci*avx_size);
                    mi[ri][ci] = m[ri][ci].markAsScalarInput();
                }
            }
        }
    }
}

template<class mmType>
void setScalarMatrixInput(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, const ScalarMatrix& m) {
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], m[ri][ci]);
        }
    }
}


template<class mmType>
void compareMatrixResults(const AVXMatrix<mmType>& mv, const ScalarMatrix& m) {
    if (mv.rows() != m.rows()) { std::cout << "m.rows,mv.rows = " << mv.rows() << " , " << m.rows() << "\n";}
    if (mv.cols() != m.cols()) { std::cout << "m.cols,mv.cols = " << mv.cols() << " , " << m.cols() << "\n";}
    for (int ri = 0; ri < mv.rows(); ++ri) {
        for (int ci = 0; ci < mv.cols(); ++ci) {
            if (abs(mmSum(mv[ri][ci])-m[ri][ci]) > 1e-10) {
                std::cout 
                    << "ri, ci, sum(mv), m = " << ri << " , " << ci << 
                    " , " << mmSum(mv[ri][ci]) << " , " << m[ri][ci] << "\n"
                ;
            }  
        }
    }
}

template<class mmType>
void compareVectorResults(const AVXVector<mmType>& vv, const ScalarVector& v) {
    if (vv.size() != v.size()) {std::cout << "m.size,mv.size = " << vv.size() << " , " << v.size() << "\n";}
    for (int ri = 0; ri < vv.size(); ++ri) {
        if (abs(mmSum(vv[ri])-v[ri]) >  1e-10) {
            std::cout << "ri,  sum(vv[ri]), v[ri] = " << ri << " , "  << mmSum(vv[ri]) << " , " << v[ri] << "\n";
        }
    }
}


template<class mmType>
void getAVXVectorGradient(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, AVXVector<mmType>& v) {
    v.resize(vi.size());
    for (int ri = 0; ri < v.size(); ++ri) {
        v[ri] = ws.diff(vi[ri]);
    }
}

template<class mmType>
void getScalarVectorGradient(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& vi, ScalarVector& v) {
    v.resize(vi.size());
    for (int ri = 0; ri < v.size(); ++ri) {
        v[ri] = ws.diff(vi[ri]);
    }
}

template<class mmType>
void getAVXMatrixGradient(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, AVXMatrix<mmType>& m) {
    m.resize(mi.rows(), mi.cols());
    for (int ri = 0; ri < m.size(); ++ri) {
        for (int ci = 0;ci < m[ri].size(); ++ci) {
            m[ri][ci] = ws.diff(mi[ri][ci]);
        }
    }
}



inline void markVectorAsInput(VectorArg& vi, const iVector& v, const bool Diff) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, v.size()));
    vi.resize(v.size());
    for (int ri = 0; ri < v.size(); ++ri) {
            v[ri].initVarAtAddress(addr_start); ++addr_start;
            if (Diff) {
                vi[ri] = v[ri].markAsInput();
            } else {
                vi[ri] = v[ri].markAsInputNoDiff();
            }   
    }
}


inline void markVectorAsOutput(VectorRes& mi, const iVector& m) {
    mi.resize(m.size());
    for (int ri = 0; ri < m.size(); ++ri) {
            mi[ri] = m[ri].markAsOutput();
    }
}

inline void markVectorAsDiff(VectorArg& vi, const iVector& v) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, v.size()));
    vi.resize(v.size());
    for (int ri = 0; ri < v.size(); ++ri) {
        v[ri].initVarAtAddress(addr_start); ++addr_start;
        vi[ri] = v[ri].markAsDiff();
    }
}

inline void markScalarVectorAsInput(ScalarVectorArg& vi, iVector& v) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, v.size()));
    vi.resize(v.size());
    for (int ri = 0; ri < v.size(); ++ri) {
        v[ri].initScalarVarAtAddress(addr_start); ++addr_start;
        vi[ri] = v[ri].markAsScalarInput();
    }
}


template<class mmType>
void setScalarVectorInput(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& vi, const ScalarVector& v) {
    for (int ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], v[ri]);
    }
}

template<class mmType, typename scalarType>  
void restructurizeData (
    const std::vector<std::vector<scalarType>>& scalar_data, 
    std::vector<AVXVector<mmType>>& mm_data
) {
    int AVX_size = aadc::mmSize<mmType>();
    mm_data.resize(scalar_data.size()/AVX_size);
    for (int i = 0; i < mm_data.size(); ++i) {
        mm_data[i].resize(scalar_data[0].size());
        for (int j=0; j<scalar_data[0].size(); j++) {
            for (int ci = 0; ci < AVX_size; ++ci) {
                toDblPtr(mm_data[i][j])[ci]=scalar_data[i*AVX_size+ci][j];
            }
        }
    }
}


template<class mmType>  
void setAVXVector(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, const AVXVector<mmType>& v) {
    for (int ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], v[ri]);
    }
}
template<class mmType>
void setScalarVectorToAVX(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, const ScalarVector& v) {
    for (int ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], mmSetConst<mmType>(v[ri]));
    }
}
template<class mmType>
void getScalarMatrixGradient(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, ScalarMatrix& m) {
    m.resize(mi.rows(), mi.cols());
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            m[ri][ci] = ws.diff(mi[ri][ci]);
        }
    }
}



template<class mmType>
void setAVXMatrixInput(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, const AVXMatrix<mmType> m) {
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], m[ri][ci]);
            //std::cout << ri << "  " << ci << " " << m[ri][ci][0] << std::endl;
        }
    }
}

template<class mmType>
void setScalarMatrixInputToAVX(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, const ScalarMatrix& m) {
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], mmSetConst<mmType>(m[ri][ci]));
        }
    }
}

} // namespace aadc
