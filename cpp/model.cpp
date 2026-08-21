#include "model.hpp"

#include "math_utils.hpp"

namespace wheelleg {

Parameters createParameters() {
    constexpr double pi = 3.14159265358979323846;
    Parameters prm{};

    // ---- base body (기준계 = base.Ai) ----
    prm.rho0p.setZero();
    prm.C00 = ang2mat(0.0, pi / 2.0, 0.0);
    prm.m0 = 111.015764646288;
    prm.J0p << 9.34382685772924, 0.0, 0.0,
               0.0, 0.27753941161572, 0.0,
               0.0, 0.0, 9.43633999493448;

    // ---- link body (기준계 = body.Ai) ----
    prm.s01p.setZero();
    prm.C01 = ang2mat(pi / 2.0, pi / 2.0, pi / 2.0);
    prm.rho1p << 6.22253616772498e-6, 0.359132273856264, 0.0;
    prm.C11 = ang2mat(pi, pi / 2.0, pi / 2.0);
    prm.m1 = 18.6342592049469;
    prm.J1p << 0.604104103452763, 1.26126031622254e-7, -3.51625272013543e-19,
               1.26126031622254e-7, 0.56572829220829, -3.57293476271758e-5,
              -3.51625272013543e-19, -3.57293476271758e-5, 0.0403172972174389;

    // ---- system ----
    prm.g = -9.80665;
    prm.F_ex = 10.0;
    prm.free_dof = {2, 7};
    return prm;
}

Eigen::Vector4d baseAttitude() {
    constexpr double pi = 3.14159265358979323846;
    return mat2ep(ang2mat(0.0, -pi / 2.0, 0.0));
}

State15 fullState(const State4& x) {
    State15 Y = State15::Zero();
    Y.segment<3>(0) << 0.0, x(0), -0.05;
    Y.segment<4>(3) = baseAttitude();
    Y(7) = x(1);
    Y.segment<3>(8) << 0.0, x(2), 0.0;
    Y.segment<3>(11).setZero();
    Y(14) = x(3);
    return Y;
}

State4 reducedState(const State15& Y) {
    State4 x;
    x << Y(1), Y(7), Y(9), Y(14);
    return x;
}

}  // namespace wheelleg
