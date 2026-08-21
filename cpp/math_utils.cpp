#include "math_utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wheelleg {

Eigen::Matrix3d tilde(const Eigen::Vector3d& x) {
    Eigen::Matrix3d out;
    out << 0.0, -x(2), x(1),
           x(2), 0.0, -x(0),
          -x(1), x(0), 0.0;
    return out;
}

Eigen::Matrix3d ang2mat(double psi, double theta, double phi) {
    Eigen::Matrix3d Rz1;
    Rz1 << std::cos(psi), -std::sin(psi), 0.0,
           std::sin(psi),  std::cos(psi), 0.0,
           0.0,            0.0,           1.0;

    Eigen::Matrix3d Rx;
    Rx << 1.0, 0.0,             0.0,
          0.0, std::cos(theta), -std::sin(theta),
          0.0, std::sin(theta),  std::cos(theta);

    Eigen::Matrix3d Rz2;
    Rz2 << std::cos(phi), -std::sin(phi), 0.0,
           std::sin(phi),  std::cos(phi), 0.0,
           0.0,            0.0,           1.0;
    return Rz1 * Rx * Rz2;
}

Eigen::Matrix3d ep2mat(const Eigen::Vector4d& p_in) {
    const Eigen::Vector4d p = p_in.normalized();
    const double e0 = p(0);
    const Eigen::Vector3d e = p.tail<3>();

    Eigen::Matrix<double, 3, 4> E;
    Eigen::Matrix<double, 3, 4> G;
    E.col(0) = -e;
    E.rightCols<3>() = tilde(e) + e0 * Eigen::Matrix3d::Identity();
    G.col(0) = -e;
    G.rightCols<3>() = -tilde(e) + e0 * Eigen::Matrix3d::Identity();
    return E * G.transpose();
}

Eigen::Vector4d mat2ep(
    const Eigen::Matrix3d& A,
    const std::optional<Eigen::Vector4d>& p_prev) {
    const double tr = A.trace();
    Eigen::Vector4d c;
    c << 1.0 + tr,
         1.0 + 2.0 * A(0, 0) - tr,
         1.0 + 2.0 * A(1, 1) - tr,
         1.0 + 2.0 * A(2, 2) - tr;

    Eigen::Index k = 0;
    const double cmax = c.maxCoeff(&k);
    if (cmax <= 0.0) {
        throw std::runtime_error("mat2ep: invalid rotation matrix");
    }
    const double s = std::sqrt(cmax);
    const double d = 0.5 / s;

    double e0 = 0.0;
    double e1 = 0.0;
    double e2 = 0.0;
    double e3 = 0.0;
    switch (k) {
        case 0:
            e0 = 0.5 * s;
            e1 = (A(2, 1) - A(1, 2)) * d;
            e2 = (A(0, 2) - A(2, 0)) * d;
            e3 = (A(1, 0) - A(0, 1)) * d;
            break;
        case 1:
            e1 = 0.5 * s;
            e0 = (A(2, 1) - A(1, 2)) * d;
            e2 = (A(0, 1) + A(1, 0)) * d;
            e3 = (A(0, 2) + A(2, 0)) * d;
            break;
        case 2:
            e2 = 0.5 * s;
            e0 = (A(0, 2) - A(2, 0)) * d;
            e1 = (A(0, 1) + A(1, 0)) * d;
            e3 = (A(1, 2) + A(2, 1)) * d;
            break;
        default:
            e3 = 0.5 * s;
            e0 = (A(1, 0) - A(0, 1)) * d;
            e1 = (A(0, 2) + A(2, 0)) * d;
            e2 = (A(1, 2) + A(2, 1)) * d;
            break;
    }

    Eigen::Vector4d p(e0, e1, e2, e3);
    p.normalize();
    if (p_prev.has_value()) {
        if (p.dot(*p_prev) < 0.0) {
            p = -p;
        }
    } else if (p(0) < 0.0) {
        p = -p;
    }
    return p;
}

double step5(double x, double x0, double h0, double x1, double h1) {
    if (x <= x0) {
        return h0;
    }
    if (x >= x1) {
        return h1;
    }
    const double a = (x - x0) / (x1 - x0);
    return h0 + (h1 - h0) * std::pow(a, 3) *
                    (10.0 - 15.0 * a + 6.0 * a * a);
}

}  // namespace wheelleg
