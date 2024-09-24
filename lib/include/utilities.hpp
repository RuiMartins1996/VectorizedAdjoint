#ifndef UTILITIES_HEADER // guard
#define UTILITIES_HEADER

// Standard Library
#include <cmath>
#include <functional>
#include <iomanip>
#include <vector>

// IO stream
#include <fstream>
#include <iostream>

// AADC library
#include <aadc/aadc.h>
#include <aadc/aadc_debug.h>

// Boost
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

const bool logs_on = true;

static void convert(std::vector<double> &alphafinal,
                    const ublas::matrix<double> &alphasmat)
{
    int Nin = alphasmat.size1();
    int Npar = alphasmat.size2();

    for (int i = 0; i < Nin; i++) {
        for (int k = 0; k < Npar; k++)
            alphafinal[Npar * i + k] = alphasmat(i, k);
    }
}

static void convert(std::vector<double> &alphafinal,
                    const std::vector<std::vector<double>> &alphasmat)
{
    int Nin = alphasmat.size();
    int Npar = alphasmat[0].size();

    for (int i = 0; i < Nin; i++) {
        for (int k = 0; k < Npar; k++)
            alphafinal[i * Npar + k] = alphasmat[i][k];
    }
}

void logavx(const aadc::mmVector<mmType> &v)
{
    int sizeOfBatch = static_cast<int>(sizeof(mmType) / sizeof(double));
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < sizeOfBatch - 1; j++)
            std::cout << v[i][j] << ",";
        std::cout << v[i][sizeOfBatch - 1] << "\n";
    }
}

template <typename State>
void log(State v)
{
    for (int i = 0; i < v.size() - 1; i++)
        std::cout << v[i] << ",";
    std::cout << v[v.size() - 1] << "\n";
};

template <typename T>
void log(ublas::matrix<T> m)
{
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2() - 1; j++) {
            std::cout << m(i, j) << "\t";
        }
        std::cout << m(i, m.size2() - 1) << "\n";
    }
}

//* Cumulative distribution function N(x)
double normalCDF(double value);

//* Functions for writing data to file
void write_vector_to_file(std::ofstream &outdata, std::vector<double> vector);
void write_vector_to_file(std::ofstream &outdata, ublas::vector<double> vector);
void write_vector_to_file(std::ofstream &outdata, ublas::matrix<double> matrix);

template <typename T>
void initialize_matrix_with_zeros(ublas::matrix<T> &m)
{
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2(); j++) {
            m(i, j) = 0.0;
        }
    }
}

double compute_euclidian_norm(ublas::vector<double> &v);

double compute_vector_sum(ublas::vector<double> &v);

double compute_vector_average(ublas::vector<double> &v);

void machineEpsilon(double EPS);

template <typename State>
double compute_error_l2(State x1, State x2)
{
    if (x1.size() != x2.size()) {
        std::cout << "Error\n";
    }

    ublas::vector<double> v1(x1.size());
    ublas::vector<double> v2(x2.size());

    for (int i = 0; i < x1.size(); i++) {
        v1(i) = x1[i];
        v2(i) = x2[i];
    }

    return ublas::norm_2(v1 - v2);
}

template <typename State>
double get_max_absolute(const State &x)
{
    std::vector<double> u(x.size());
    for (int i = 0; i < x.size(); i++) {
        u[i] = std::abs(x[i]);
    }

    double max = 0;
    for (auto &e : u) {
        max = std::max(max, e);
    }

    return max;
}

template <typename State>
double compute_error(State x1, State x2)
{
    if (x1.size() != x2.size()) {
        std::cout << "Error\n";
    }

    ublas::vector<double> v1(x1.size());
    ublas::vector<double> v2(x2.size());

    for (int i = 0; i < x1.size(); i++) {
        v1(i) = x1[i];
        v2(i) = x2[i];
    }

    return ublas::norm_inf(v1 - v2);
}

void get_number_error_pair(double &d, double &e, int &n)
{
    if (e >= 1) {
        int m = (int)std::ceil(
            std::log10(e)); // this is the position of the
                            // first significant digit before the decimal point
        e = e * std::pow(10, n - m);
        e = std::round(e);
        e = e * std::pow(10, m - n);

        d = d * std::pow(10, n - m);
        d = std::round(d);
        d = d * std::pow(10, m - n);

        if (n <= m) {
            n = 0;
        } else {
            n = n - m;
        }

    } else {
        int m = (int)std::ceil(
            -std::log10(e)); // this is the position of the first significant
                             // digit after the decimal point
        e = e * std::pow(10, m + n - 1);
        e = std::round(e);
        e = e / std::pow(10, m + n - 1);

        d = d * std::pow(10, m + n - 1);
        d = std::round(d);
        d = d / std::pow(10, m + n - 1);

        n = m + n - 1;
        ;
    }
}

void write_number_error(double d, double e, int n)
{
    get_number_error_pair(d, e, n);

    std::cout << std::fixed << std::setprecision(n) << d << ","
              << std::setprecision(n) << e;
}

void write_number_error_pair(double d, double e, int n, std::ofstream &file)
{
    get_number_error_pair(d, e, n);

    file << std::fixed << std::setprecision(n) << d << ","
         << std::setprecision(n) << e;
}

int first_significant(double d)
{
    if (d >= 1) {
        double intpart;
        double fractpart = modf(d, &intpart);
        return (int)std::ceil(-std::log10(fractpart));
    } else {
        return (int)std::ceil(-std::log10(d));
    }
}

double compute_std_deviation(std::vector<double> v, double v_mean)
{
    int N = v.size();

    double sum = 0.0;

    for (int i = 0; i < N; i++) {
        sum += (v[i] - v_mean) * (v[i] - v_mean);
    }

    if (N > 1) {
        return std::sqrt(sum / static_cast<double>(N - 1));
    } else {
        return std::sqrt(sum / static_cast<double>(N));
    }
}

#endif