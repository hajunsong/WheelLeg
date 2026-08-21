#pragma once

#include <Eigen/Dense>

#include <array>

namespace wheelleg {

using State15 = Eigen::Matrix<double, 15, 1>;
using State4 = Eigen::Matrix<double, 4, 1>;
using Matrix7d = Eigen::Matrix<double, 7, 7>;
using Vector7d = Eigen::Matrix<double, 7, 1>;
using Matrix4d = Eigen::Matrix4d;
using RowVector4d = Eigen::Matrix<double, 1, 4>;

struct Parameters {
    Eigen::Vector3d rho0p;
    Eigen::Matrix3d C00;
    double m0;
    Eigen::Matrix3d J0p;

    Eigen::Vector3d s01p;
    Eigen::Matrix3d C01;
    Eigen::Vector3d rho1p;
    Eigen::Matrix3d C11;
    double m1;
    Eigen::Matrix3d J1p;

    double g;
    double F_ex;
    std::array<int, 2> free_dof;  // MATLAB과 같은 1-base 번호
};

Parameters createParameters();
Eigen::Vector4d baseAttitude();
State15 fullState(const State4& x);
State4 reducedState(const State15& Y);

}  // namespace wheelleg
