#pragma once

#include <Eigen/Dense>

#include <optional>

namespace wheelleg {

Eigen::Matrix3d tilde(const Eigen::Vector3d& x);
Eigen::Matrix3d ang2mat(double psi, double theta, double phi);
Eigen::Matrix3d ep2mat(const Eigen::Vector4d& p);
Eigen::Vector4d mat2ep(
    const Eigen::Matrix3d& A,
    const std::optional<Eigen::Vector4d>& p_prev = std::nullopt);
double step5(double x, double x0, double h0, double x1, double h1);

}  // namespace wheelleg
