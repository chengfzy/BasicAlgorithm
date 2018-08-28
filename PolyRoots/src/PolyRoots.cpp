#include "PolyRoots.h"

// Compute the roots based on polynominal coefficients
void PolyRoots::compute(const Eigen::VectorXd& coeffs) {
    // find the non-zero range [n0, n1]
    Eigen::Index kN = coeffs.size();
    Eigen::Index n0{0}, n1{coeffs.size() - 1};
    while (n0 < kN && coeffs[n0] == 0) {
        ++n0;
    }
    while (n1 >= n0 && coeffs[n1] == 0) {
        --n1;
    }

    // initialize roots
    roots_.resize(kN - n0 - 1);
    roots_.setZero();

    // construct matrix A and compute eigenvalues
    const Eigen::Index kNonZeroN = n1 - n0;
    if (kNonZeroN > 0) {
        Eigen::MatrixXd A(kNonZeroN, kNonZeroN);
        A.setZero();
        A.diagonal(-1).setOnes();
        A.row(0) = -coeffs.segment(n0 + 1, kNonZeroN) / coeffs[n0];
        solver_.compute(A);
        roots_.tail(kNonZeroN) = solver_.eigenvalues();
    }
}

// Return the real roots
std::vector<double> PolyRoots::realRoots(const double& absImagThresh) const {
    std::vector<double> roots;
    for (Eigen::Index i = 0; i < roots_.size(); ++i) {
        if (std::abs(roots_[i].imag()) < absImagThresh) {
            roots.emplace_back(roots_[i].real());
        }
    }
    return roots;
}