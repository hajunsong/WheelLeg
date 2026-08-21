#pragma once

#include "model.hpp"

#include <Eigen/Dense>

#include <functional>

namespace wheelleg {

struct LinearModel {
    Matrix4d A;
    Eigen::Matrix<double, 4, 1> B;
};

struct LqrResult {
    Matrix4d P;
    RowVector4d K;
    double care_relative_residual;
    double u1_condition_number;
    double imaginary_ratio;
    double symmetry_ratio;
};

State4 nonlinearDynamics(
    double t,
    const State4& x,
    const Parameters& prm,
    double force);
State4 uprightEquilibrium(const Parameters& prm, double y_ref = 0.0);
LinearModel numericalLinearization(
    const std::function<State4(const State4&, double)>& dynamics,
    const State4& xe,
    double ue);
Eigen::Matrix4d controllabilityMatrix(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B);
bool isStabilizable(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B,
    double tolerance = 1e-9);
LqrResult solveCareHamiltonian(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B,
    const Matrix4d& Q,
    double R);
double angleError(double q, double q_ref);

}  // namespace wheelleg
