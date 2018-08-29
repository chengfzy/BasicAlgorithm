#pragma once
#include <algorithm>
#include <array>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SparseLU"

/**
 * @brief Cubic spline to fit data. For input data (t, X), which t is Nx1 parameter data, x is N x (Dim - 1) data, try
 * to find cubic spine fi(t) = ai * (t - ti)^3 + bi * (t - ti)^2 + ci * (t - ti) + xi, where i in [0, N-1]. If not input
 * parameter data "t", will calculate chord length of x as parameter.
 * @tparam Dim  Data dimension
 */
template <unsigned int Dim = 2>
class CubicSpline {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
   public:
    /**
     * @brief Constructor
     */
    CubicSpline();

   public:
    /**
     * @brief   Build cubic spline from input data
     * @tparam N    Data length
     * @param data  Input data which include parameter vector "t"
     */
    template <int N = Eigen::Dynamic>
    void build(const Eigen::Matrix<double, N, Dim>& data);

    /**
     * @brief Build cubic spline from input data without parameter column "t", could not use same function name "build"
     * because not suitable for MatrixXd input
     * @tparam N    Data length
     * @param data  Input data which not include parameter vector "t"
     */
    template <int N = Eigen::Dynamic>
    void build(const Eigen::Matrix<double, N, Dim - 1>& data);

    /**
     * @brief Evaluate value at given position t
     * @param t     Given position
     * @return Spline value f(t) at position t
     */
    Eigen::Matrix<double, Dim - 1, 1> operator()(const double& t) const;

    /**
     * @brief Evaluate derive f'(t) at given position t
     * @param t Given position
     * @return Derive value f'(t) at position
     */
    Eigen::Matrix<double, Dim - 1, 1> derive(const double& t) const;

    /**
     * @brief Get coefficients at index
     * @param index Index, should less than N
     * @return Coefficients, i,e, [a, b, c, d]
     */
    Eigen::Matrix<double, 4, Dim - 1> coeffAt(const std::size_t& index) const;

    /**
     * @brief Return the start parameters ti at index
     * @param index Index of the spline
     * @return Start parameters ti
     */
    inline const double& startParamsAt(const std::size_t& index) const { return data_(index, 0); }

   private:
    /**
     * @brief Normalize parameter vector
     */
    inline double normalizeParamVector(const double& t) const { return (t - tMin_) / (tMax_ - tMin_); }

    /**
     * @brief Compute parameter vector "t" if input data not included
     * @tparam N    Data length
     * @param data  Input data which not include parameter vector "t"
     */
    template <int N = Eigen::Dynamic>
    void calParamVector(const Eigen::Matrix<double, N, Dim - 1>& data);

    /**
     * @brief Build cubic spline
     */
    void build();

   private:
    Eigen::MatrixXd data_;  // data include parameter vector, (t, x)
    double tMin_;           // min value of parameter vector
    double tMax_;           // max value of parameter vector
    // poly coefficients, a, b, c, and d = y, the last coefficient a[n-1], b[n-1], c[n-1] for right extrapolation
    std::array<Eigen::MatrixX3d, Dim - 1> coeffs_;
};

// Type Definition
using CubicSpline1d = CubicSpline<1>;
using CubicSpline2d = CubicSpline<2>;
using CubicSpline3d = CubicSpline<3>;

// Constructor
template <unsigned int Dim>
CubicSpline<Dim>::CubicSpline() : tMin_(0.0), tMax_(1.0) {
    static_assert(Dim > 1, "Should at Least Have 2 Dimension Data for Cubic Spline, and Forbid to Use Dynamic Size");
}

// Build cubic spline from input data
template <unsigned int Dim>
template <int N>
void CubicSpline<Dim>::build(const Eigen::Matrix<double, N, Dim>& data) {
    // check data length
    const Eigen::Index kN = data.rows();
    if (kN < 3) {
        throw std::length_error("Data Length Should More than 3 for Cubic Spline");
    }

    // check parameter column t, t0 < t1 < ... < tn
    for (Eigen::Index i = 0; i < kN - 1; ++i) {
        assert(data(i, 0) < data(i + 1, 0));
    }

    // initialize data
    data_ = data;

    // normalize parameter vector t
    tMin_ = data(0, 0);
    tMax_ = data(kN - 1, 0);
    for (Eigen::Index i = 0; i < kN; ++i) {
        data_(i, 0) = normalizeParamVector(data_(i, 0));
    }

    // build spline
    build();
}

// Build cubic spline from input data without parameter column "t"
template <unsigned int Dim>
template <int N>
void CubicSpline<Dim>::build(const Eigen::Matrix<double, N, Dim - 1>& data) {
    // initialize data: allocate memory
    data_.resize(data.rows(), Dim);

    // calculate parameter vector and assign data
    calParamVector(data);
    data_.rightCols(Dim - 1) = data;
    tMin_ = 0.0;
    tMax_ = 1.0;

    // build spline
    build();
}

// Evaluate value at given position f(t)
template <unsigned int Dim>
Eigen::Matrix<double, Dim - 1, 1> CubicSpline<Dim>::operator()(const double& t) const {
    const Eigen::Index kN = data_.rows();

    // normalization
    double tt = normalizeParamVector(t);

    // find the index t[idx] < t, and idx = 0 if t < t0; idx = n - 1 if t > t[n-1]
    Eigen::Index idx{0};
    while (idx < kN && data_(idx, 0) < tt) {
        ++idx;
    }
    if (idx > 0) --idx;

    double h = tt - data_(idx, 0);
    Eigen::Matrix<double, Dim - 1, 1> ret = Eigen::Matrix<double, Dim - 1, 1>::Zero();
    if (tt < 0) {
        // left extrapolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = (coeffs_[i](0, 1) * h + coeffs_[i](0, 2)) * h + data_(0, i + 1);
        }
    } else if (tt > 1) {
        // right extrapolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = (coeffs_[i](kN - 1, 1) * h + coeffs_[i](kN - 1, 2)) * h + data_(kN - 1, i + 1);
        }
    } else {
        // interpolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = ((coeffs_[i](idx, 0) * h + coeffs_[i](idx, 1)) * h + coeffs_[i](idx, 2)) * h + data_(idx, i + 1);
        }
    }
    return ret;
}

// Evaluate derive f'(t) at given position t
template <unsigned int Dim>
Eigen::Matrix<double, Dim - 1, 1> CubicSpline<Dim>::derive(const double& t) const {
    const Eigen::Index kN = data_.rows();

    // normalization
    double tt = normalizeParamVector(t);

    // find the index t[idx] < t, and idx = 0 if t < t0; idx = n - 1 if t > t[n-1]
    Eigen::Index idx{0};
    while (idx < kN && data_(idx, 0) < tt) {
        ++idx;
    }
    if (idx > 0) --idx;

    double h = tt - data_(idx, 0);
    Eigen::Matrix<double, Dim - 1, 1> ret = Eigen::Matrix<double, Dim - 1, 1>::Zero();
    if (tt < 0) {
        // left extrapolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = 2 * coeffs_[i](0, 1) * h + coeffs_[i](0, 2);
        }
    } else if (tt > 1) {
        // right extrapolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = 2 * coeffs_[i](kN - 1, 1) * h + coeffs_[i](kN - 1, 2);
        }
    } else {
        // interpolation
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            ret[i] = (3 * coeffs_[i](idx, 0) * h + 2 * coeffs_[i](idx, 1)) * h + coeffs_[i](idx, 2);
        }
    }
    return ret / (tMax_ - tMin_);
}

// Get coefficients at index
template <unsigned int Dim>
Eigen::Matrix<double, 4, Dim - 1> CubicSpline<Dim>::coeffAt(const std::size_t& index) const {
    Eigen::Matrix<double, 4, Dim - 1> coeff = Eigen::Matrix<double, 4, Dim - 1>::Zero();
    if (index < data_.rows()) {
        for (std::size_t i = 0; i < Dim - 1; ++i) {
            coeff.col(i).topRows(3) = coeffs_[i].row(index);  // a, b, c
            coeff(3, i) = data_(index, i + 1);                // d = x
        }
    }
    return coeff;
}

// Compute parameter vector "t" if input data not included
template <unsigned int Dim>
template <int N>
void CubicSpline<Dim>::calParamVector(const Eigen::Matrix<double, N, Dim - 1>& data) {
    const Eigen::Index kN = data.rows();

    // compute chord length and partial sum
    data_.block(1, 0, kN - 1, 1) = (data.topRows(kN - 1) - data.bottomRows(kN - 1)).rowwise().norm();
    std::partial_sum(data_.data(), data_.data() + kN, data_.data());

    // normalization
    data_.col(0) /= data_(kN - 1, 0);
    data_(kN - 1, 0) = 1.0;  // force assign to 1
}

// Build cubic spline
template <unsigned int Dim>
void CubicSpline<Dim>::build() {
    using Triplet = Eigen::Triplet<double>;
    const Eigen::Index kN = data_.rows();

    // initialize coefficients
    for (auto& c : coeffs_) {
        c.resize(kN, 3);
        c.setZero();
    }

    // extract t and calculate h and dx/h, h[i] = t[i+1] - t[i-1], dxh = (x[i+1] - x[i]) / h[i]
    const auto& t = data_.col(0);
    auto h = t.bottomRows(kN - 1) - t.topRows(kN - 1);
    auto dxhMat =
        (data_.block(1, 1, kN - 1, Dim - 1) - data_.block(0, 1, kN - 1, Dim - 1)).array().colwise() / h.array();

    // calculate b for each dimension data, A * b = a
    for (Eigen::Index n = 0; n < Dim - 1; ++n) {
        const auto& dxh = dxhMat.col(n);
        Eigen::MatrixX3d& coeff = coeffs_[n];

        // set up matrix A
        Eigen::SparseMatrix<double> A(kN, kN);
        Eigen::Matrix<double, Eigen::Dynamic, 1> a(kN, 1);
        std::vector<Eigen::Triplet<double>> triplets;
        for (Eigen::Index i = 1; i < kN - 1; ++i) {
            triplets.emplace_back(Triplet(i, i - 1, h[i - 1] / 3.0));
            triplets.emplace_back(Triplet(i, i, 2.0 / 3.0 * (h[i - 1] + h[i])));
            triplets.emplace_back(i, i + 1, h[i] / 3.0);
            a[i] = dxh[i] - dxh[i - 1];
        }

        // boundary conditions, for extrapolation we use second order polynomials, i.e.,
        // f(t) = b0 * (t - t0)^2 + c0 * (t - x0) + x0, for t <= t0
        // f(t) = b[n-1] * (t - t[n-1])^2 + c[n-1] * (t - t[n-1] + x[n-1], for t >= t[n-1]
        // so f'' = 2 * b = 0
        triplets.emplace_back(Triplet(0, 0, 2.0));
        triplets.emplace_back(Triplet(kN - 1, kN - 1, 2.0));
        a[0] = 0;
        a[kN - 1] = 0;

        // set data to A
        A.setFromTriplets(triplets.begin(), triplets.end());

        // solve b from equation A * b = a
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        coeff.col(1) = solver.solve(a);

        // solve a and c
        for (Eigen::Index i = 0; i < kN - 1; ++i) {
            coeff(i, 0) = (coeff(i + 1, 1) - coeff(i, 1)) / h[i] / 3.0;
            coeff(i, 2) = dxh[i] - (coeff(i + 1, 1) + 2.0 * coeff(i, 1)) * h[i] / 3.0;
        }
        // right boundary: fn(t) = bn * (t - tn)^2 + cn * (t - tn) + tn, for t >= tn
        // fn'(tn) = f[n-1]'(tn), so cn = 3 * a[n-1] * (tn - t[n-1])^2 + 2 * b[n-1] * (tn - t[n-1]) + c[n-1]
        coeff(kN - 1, 2) =
            3.0 * coeff(kN - 2, 0) * h[kN - 2] * h[kN - 2] + 2.0 * coeff(kN - 2, 1) * h[kN - 2] + coeff(kN - 2, 2);
    }
}
