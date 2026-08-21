#include "dYdt.hpp"

#include "math_utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wheelleg {

State15 dYdt(
    double t,
    const State15& Y,
    const Parameters& prm,
    std::optional<double> control_force,
    DyDtOut* out) {
    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    using Vector6d = Eigen::Matrix<double, 6, 1>;

    // ---------------- Y2qdq ----------------
    const Eigen::Vector3d r0 = Y.segment<3>(0);
    Eigen::Vector4d p0 = Y.segment<4>(3);
    const double q1 = Y(7);
    const Eigen::Vector3d dr0 = Y.segment<3>(8);
    const Eigen::Vector3d w0 = Y.segment<3>(11);
    const double dq1 = Y(14);

    p0.normalize();  // RK4 중간 stage를 단위 quaternion에 투영

    // ---------------- base body ----------------
    const double e0 = p0(0);
    const Eigen::Vector3d e = p0.tail<3>();
    Eigen::Matrix<double, 3, 4> E0;
    Eigen::Matrix<double, 3, 4> G0;
    E0.col(0) = -e;
    E0.rightCols<3>() = tilde(e) + e0 * Eigen::Matrix3d::Identity();
    G0.col(0) = -e;
    G0.rightCols<3>() = -tilde(e) + e0 * Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d A0 = E0 * G0.transpose();

    const Eigen::Matrix3d AC00 = A0 * prm.C00;
    const Eigen::Matrix3d J0c = AC00 * prm.J0p * AC00.transpose();

    const Eigen::Vector3d rho0 = A0 * prm.rho0p;
    const Eigen::Vector3d r0c = r0 + rho0;

    const Eigen::Matrix3d w0t = tilde(w0);
    const Eigen::Matrix3d r0t = tilde(r0);
    const Eigen::Vector3d dr0c = dr0 + w0t * rho0;

    const Eigen::Matrix3d dr0t = tilde(dr0);
    const Eigen::Matrix3d dr0ct = tilde(dr0c);
    const Eigen::Matrix3d r0ct = tilde(r0c);

    Vector6d Y0h;
    Y0h << dr0 + r0t * w0, w0;

    // 외력: RecurDyn TRANSLATIONAL_FORCE, global FY [N]
    const double Fy = control_force.has_value()
                          ? *control_force
                          : step5(t, 0.0, prm.F_ex, 1.0, -prm.F_ex);
    const Eigen::Vector3d f0c(0.0, Fy, prm.m0 * prm.g);
    const Eigen::Vector3d t0c = Eigen::Vector3d::Zero();

    Matrix6d M0h;
    M0h.topLeftCorner<3, 3>() = prm.m0 * Eigen::Matrix3d::Identity();
    M0h.topRightCorner<3, 3>() = -prm.m0 * r0ct;
    M0h.bottomLeftCorner<3, 3>() = prm.m0 * r0ct;
    M0h.bottomRightCorner<3, 3>() = J0c - prm.m0 * r0ct * r0ct;

    Vector6d Q0h;
    Q0h << f0c + prm.m0 * dr0ct * w0,
           t0c + r0ct * f0c + prm.m0 * r0ct * dr0ct * w0 -
               w0t * J0c * w0;

    // ---------------- link body ----------------
    Eigen::Matrix3d A01pp;
    A01pp << std::cos(q1), -std::sin(q1), 0.0,
             std::sin(q1),  std::cos(q1), 0.0,
             0.0,           0.0,          1.0;
    const Eigen::Matrix3d A1 = A0 * prm.C01 * A01pp;
    const Eigen::Vector3d s01 = A1 * prm.s01p;
    const Eigen::Vector3d r1 = r0 + s01;

    const Eigen::Matrix3d AC11 = A1 * prm.C11;
    const Eigen::Matrix3d J1c = AC11 * prm.J1p * AC11.transpose();

    const Eigen::Vector3d rho1 = A1 * prm.rho1p;
    const Eigen::Vector3d r1c = r1 + rho1;
    const Eigen::Matrix3d r1ct = tilde(r1c);

    const Eigen::Vector3d H1 = A0 * prm.C01 * Eigen::Vector3d::UnitZ();
    const Eigen::Vector3d w1 = w0 + H1 * dq1;
    const Eigen::Matrix3d w1t = tilde(w1);
    const Eigen::Vector3d dr1 = dr0 + w0t * s01;
    const Eigen::Matrix3d r1t = tilde(r1);

    Vector6d B1;
    B1 << r1t * H1, H1;

    const Eigen::Matrix3d dr1t = tilde(dr1);
    const Eigen::Vector3d dr1c = dr1 + w1t * rho1;
    const Eigen::Matrix3d dr1ct = tilde(dr1c);

    const Eigen::Vector3d dH1 = w0t * H1;
    Vector6d D1;
    D1 << dr1t * H1 + r1t * dH1, dH1;
    D1 *= dq1;
    const Vector6d Y1h = Y0h + B1 * dq1;

    const Eigen::Vector3d f1c(0.0, 0.0, prm.m1 * prm.g);
    const Eigen::Vector3d t1c = Eigen::Vector3d::Zero();

    Matrix6d M1h;
    M1h.topLeftCorner<3, 3>() = prm.m1 * Eigen::Matrix3d::Identity();
    M1h.topRightCorner<3, 3>() = -prm.m1 * r1ct;
    M1h.bottomLeftCorner<3, 3>() = prm.m1 * r1ct;
    M1h.bottomRightCorner<3, 3>() = J1c - prm.m1 * r1ct * r1ct;

    Vector6d Q1h;
    Q1h << f1c + prm.m1 * dr1ct * w1,
           t1c + r1ct * f1c + prm.m1 * r1ct * dr1ct * w1 -
               w1t * J1c * w1;

    // ---------------- mass / force ----------------
    const Matrix6d K1 = M1h;
    const Matrix6d K0 = K1 + M0h;
    const Vector6d L1 = Q1h;
    const Vector6d L0 = L1 + Q0h - K1 * D1;

    // ---------------- EQM ----------------
    Matrix7d M;
    M.topLeftCorner<6, 6>() = K0;
    M.topRightCorner<6, 1>() = K1 * B1;
    M.bottomLeftCorner<1, 6>() = B1.transpose() * K1;
    M(6, 6) = (B1.transpose() * K1 * B1)(0, 0);

    Vector7d Q;
    Q.head<6>() = L0;
    Q(6) = (B1.transpose() * (L1 - K1 * D1))(0, 0);

    for (int i = 0; i < 7; ++i) {
        const int matlab_dof = i + 1;
        const bool free = std::find(
                              prm.free_dof.begin(),
                              prm.free_dof.end(),
                              matlab_dof) != prm.free_dof.end();
        if (!free) {
            M.row(i).setZero();
            M(i, i) = 1.0;
            Q(i) = 0.0;
        }
    }

    const Eigen::FullPivLU<Matrix7d> lu(M);
    if (!lu.isInvertible()) {
        throw std::runtime_error("dYdt: singular constrained mass matrix");
    }
    const Vector7d ddq = lu.solve(Q);
    const Vector6d dY0h = ddq.head<6>();
    const double ddq1 = ddq(6);

    // ---------------- base body acceleration ----------------
    const Eigen::Vector4d dp0 = 0.5 * E0.transpose() * w0;
    Matrix6d T0 = Matrix6d::Zero();
    T0.topLeftCorner<3, 3>().setIdentity();
    T0.topRightCorner<3, 3>() = -r0t;
    T0.bottomRightCorner<3, 3>().setIdentity();
    Vector6d R0 = Vector6d::Zero();
    R0.head<3>() = dr0t * w0;
    const Vector6d dY0b = T0 * dY0h - R0;
    const Eigen::Vector3d ddr0 = dY0b.head<3>();
    const Eigen::Vector3d dw0 = dY0b.tail<3>();

    const Eigen::Matrix3d dw0t = tilde(dw0);
    const Eigen::Vector3d ddr0c =
        ddr0 + dw0t * rho0 + w0t * w0t * rho0;

    // ---------------- link body acceleration ----------------
    const Vector6d dY1h = dY0h + B1 * ddq1 + D1;
    Matrix6d T1 = Matrix6d::Zero();
    T1.topLeftCorner<3, 3>().setIdentity();
    T1.topRightCorner<3, 3>() = -r1t;
    T1.bottomRightCorner<3, 3>().setIdentity();
    Matrix6d dT1 = Matrix6d::Zero();
    dT1.topRightCorner<3, 3>() = -dr1t;
    const Vector6d dY1b = dT1 * Y1h + T1 * dY1h;
    const Eigen::Vector3d ddr1 = dY1b.head<3>();
    const Eigen::Vector3d dw1 = dY1b.tail<3>();
    const Eigen::Matrix3d dw1t = tilde(dw1);
    const Eigen::Vector3d ddr1c =
        ddr1 + dw1t * rho1 + w1t * w1t * rho1;

    // ---------------- dqddq2Yp ----------------
    State15 Yp;
    Yp << dr0, dp0, dq1, ddr0, dw0, ddq1;

    if (out != nullptr) {
        out->A0 = A0;
        out->A1 = A1;
        out->r0c = r0c;
        out->r1c = r1c;
        out->dr0c = dr0c;
        out->dr1c = dr1c;
        out->ddr0c = ddr0c;
        out->ddr1c = ddr1c;
        out->w1 = w1;
        out->dw0 = dw0;
        out->dw1 = dw1;
        out->ddq1 = ddq1;
        out->Fy = Fy;
        out->M = M;
        out->Q = Q;
    }
    return Yp;
}

}  // namespace wheelleg
