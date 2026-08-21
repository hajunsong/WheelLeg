#pragma once

#include "model.hpp"

#include <optional>

namespace wheelleg {

struct DyDtOut {
    Eigen::Matrix3d A0;
    Eigen::Matrix3d A1;
    Eigen::Vector3d r0c;
    Eigen::Vector3d r1c;
    Eigen::Vector3d dr0c;
    Eigen::Vector3d dr1c;
    Eigen::Vector3d ddr0c;
    Eigen::Vector3d ddr1c;
    Eigen::Vector3d w1;
    Eigen::Vector3d dw0;
    Eigen::Vector3d dw1;
    double ddq1;
    double Fy;
    Matrix7d M;
    Vector7d Q;
};

State15 dYdt(
    double t,
    const State15& Y,
    const Parameters& prm,
    std::optional<double> control_force = std::nullopt,
    DyDtOut* out = nullptr);

}  // namespace wheelleg
