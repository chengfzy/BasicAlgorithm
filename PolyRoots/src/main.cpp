/**
 * Calculate the polynominal roots. f(x) = an * x^n + a(n-1) * x^(n-1) + ... + a1 * x + a0
 */

#include <iostream>
#include "PolyRoots.h"

using namespace std;
using namespace Eigen;

// print polynominal expression
void printPoly(const Eigen::VectorXd& coeff) {
    const std::size_t N = coeff.size();
    cout << "f(x) = ";
    for (std::size_t i = 0; i < N; ++i) {
        // skill 0 * x^n
        if (coeff[i] == 0) continue;

        // print a[i]
        if (0 == i) {
            cout << coeff[i];
        } else {
            cout << abs(coeff[i]);
        }

        // print x^(N-1-i)
        int deg = N - 1 - i;
        if (deg > 1) {
            cout << " * x^" << deg;
        } else if (deg == 1) {
            cout << " * x";
        }

        // print +/-
        if (N - 1 != i) {
            cout << (coeff[i + 1] > 0 ? " + " : " - ");
        }
    }
    cout << endl;
}

// operator << overload for Eigen::VectorXcd
std::ostream& operator<<(std::ostream& out, const Eigen::VectorXcd& vec) {
    for (std::size_t i = 0; i < vec.size(); ++i) {
        out << vec[i];
        if (i != vec.size() - 1) {
            out << ", ";
        }
    }
    return out;
}

// operator << overload for vector<T>
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    for (std::size_t i = 0; i < vec.size(); ++i) {
        out << vec[i];
        if (i != vec.size() - 1) {
            out << ", ";
        }
    }
    return out;
}

int main(int argc, char* argv[]) {
    PolyRoots polyRoots;

    // calculate the root of f(x) = 3x^2 - 2x - 4;
    Vector3d c1(3, -2, -4);
    polyRoots.compute(c1);
    printPoly(c1);
    cout << "roots = " << polyRoots.roots() << endl;
    cout << "real roots = " << polyRoots.realRoots() << endl << endl;

    // calculate the root of f(x) = 3x^3 - 2x^2 - 4x + 5;
    Vector4d c2(3, -2, -4, 5);
    polyRoots.compute(c2);
    printPoly(c2);
    cout << "roots = " << polyRoots.roots() << endl;
    cout << "real roots = " << polyRoots.realRoots() << endl << endl;

    return 0;
}