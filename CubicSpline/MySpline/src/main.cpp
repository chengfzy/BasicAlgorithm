#include <fstream>
#include <iomanip>
#include <iostream>
#include "CubicSpline.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"
#include "gflags/gflags.h"
#include "glog/logging.h"

using namespace std;
using namespace Eigen;
using namespace boost;
using namespace boost::algorithm;

// Tokenize string with specified separator
std::vector<string> tokenize(const std::string& line, const std::string& separator = ",") {
    using Tokenizer = tokenizer<char_separator<char>>;
    char_separator<char> sep{separator.c_str()};
    Tokenizer tok(line, sep);
    vector<string> tokens;
    tokens.assign(tok.begin(), tok.end());
    return tokens;
}

// simple test, similar as tkSpline/main.cpp, the only different is the normlization of t
void simpleTest() {
    MatrixX2d data(5, 2);
    data << 0.100000, 0.100000, 0.400000, 0.700000, 1.200000, 0.600000, 1.800000, 1.100000, 2.000000, 0.900000;
    cout << "data = " << endl << data << endl << endl;

    // build and then interpolation
    CubicSpline2d spline;
    spline.build(data);

    for (int i = -50; i < 260; i += 50) {
        double t = 0.01 * i;
        cout << setprecision(15) << t << ", " << spline(t).transpose() << ", " << spline.derive(t) << endl;
    }
}

// generate data of sin(t), and fit it then interpolate and evaluate derive
void sinTest() {
    auto func = [](const double& t) { return sin(t); };
    auto dfunc = [](const double& t) { return cos(t); };

    // generate data f(t) = sin(t), t = [-1.0, 2.0]
    const double t0 = -1.0;
    const double t1 = 2.0;
    const double dt = 0.05;
    const size_t kN = static_cast<size_t>((t1 - t0) / dt) + 1;
    MatrixX2d data(kN, 2);
    for (size_t i = 0; i < kN; ++i) {
        double t = t0 + i * dt;
        double x = func(t);
        data(i, 0) = t;
        data(i, 1) = x;
    }

    // spine fit
    CubicSpline2d spline;
    spline.build(data);

    // interpolation and evaluate derive
    for (int i = -50; i < 251; i += 10) {
        double t = 0.01 * i;
        double x = func(t);
        double dx = dfunc(t);
        auto fx = spline(t);
        auto dfx = spline.derive(t);
        cout << setprecision(15) << "t = " << t << ", x = " << x << ", fx = " << fx.transpose() << ", dx = " << dx
             << ", dfx = " << dfx.transpose() << endl;
    }
}

// complex test: used to fit GPS data
void complexTest() {
    // load center line local data
    const int kDataSize = 700;
    MatrixX3d centerLineData(kDataSize, 3);
    fstream inFile("../../data/CenterLineLocal.txt", ios::in);
    if (!inFile.is_open()) {
        cout << "cannot open input file." << endl;
    }
    // read data: longitude, latitude, altitude
    string line;
    int lineCount{0};
    while (getline(inFile, line) && lineCount < kDataSize) {
        vector<string> tokens = tokenize(line, ",");
        if (tokens.size() != 3) {
            break;
        }
        centerLineData(lineCount, 0) = lexical_cast<double>(trim_copy(tokens[0]));
        centerLineData(lineCount, 1) = lexical_cast<double>(trim_copy(tokens[1]));
        centerLineData(lineCount, 2) = lexical_cast<double>(trim_copy(tokens[2]));
        ++lineCount;
    }

    // only fit x, y data
    MatrixX2d data = centerLineData.leftCols(2);
    cout << "data size = " << data.rows() << " x " << data.cols() << endl;
    CubicSpline3d spline;
    spline.build(data);

    // interpolation to evaluate data and save to file
    fstream file("spline.txt", ios::out);
    if (!file.is_open()) {
        cout << "cannot open or create file" << endl;
    }
    const double kStart{-0.05};
    const double kEnd{1.05};
    const double kStep{0.001};
    const int kInterpSize = static_cast<int>(ceil((kEnd - kStart) / kStep));
    for (size_t i = 0; i < kInterpSize; ++i) {
        double t = kStart + kStep * i;
        Vector2d interpPos = spline(t);
        file << setprecision(15) << t << ", " << interpPos[0] << ", " << interpPos[1] << endl;
    }
    file.close();
}

int main(int argc, char* argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::ParseCommandLineFlags(&argc, &argv, true);

    //    simpleTest();
    sinTest();
    //     complexTest();

    google::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();
    return 0;
}