#pragma once
#include "Eigen/Core"
#include "Eigen/Eigenvalues"

/**
 * @brief Compute the Roots of Polynominal. f(x) = a[n] * x^n + a[n-1] * x^(n-1) + ... + a1 * x + a0
 */
class PolyRoots {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
   public:
    PolyRoots() = default;

    /**
     * @brief Compute the roots based on polynominal coefficients
     * @param coeffs    Polynominal coefficients [a(n), a(n-1), ..., a1, a0]
     */
    void compute(const Eigen::VectorXd& coeffs);

    /**
     * @brief Return the roots
     * @return Roots
     */
    inline const Eigen::VectorXcd& roots() const { return roots_; }

    /**
     * @brief Return the real roots
     * @param absImagThresh the maximum bound of the imaginary part of a complex
     * @return Real roots
     */
    std::vector<double> realRoots(const double& absImagThresh = Eigen::NumTraits<double>::dummy_precision()) const;

   private:
    Eigen::EigenSolver<Eigen::MatrixXd> solver_;  // eigen solver
    Eigen::VectorXcd roots_;                      // roots
};