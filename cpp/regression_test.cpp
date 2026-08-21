#include "lqr.hpp"
#include "model.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

int main() {
    using namespace wheelleg;

    const Parameters prm = createParameters();
    const State4 xe = uprightEquilibrium(prm);
    const auto dynamics = [&prm](const State4& x, double force) {
        return nonlinearDynamics(0.0, x, prm, force);
    };
    if (dynamics(xe, 0.0).norm() > 1e-10) {
        std::cerr << "equilibrium residual regression failed\n";
        return 1;
    }

    const LinearModel linear = numericalLinearization(dynamics, xe, 0.0);
    if (std::abs(linear.A(2, 1) - 1.272523158) > 1e-8 ||
        std::abs(linear.A(3, 1) - 24.65311081) > 1e-7 ||
        std::abs(linear.B(2) - 0.008713930133) > 1e-10 ||
        std::abs(linear.B(3) - 0.01939002867) > 1e-10) {
        std::cerr << "A/B regression failed\n";
        return 1;
    }

    constexpr double pi = 3.14159265358979323846;
    const Eigen::Vector4d maximum_state(0.2, 5.0 * pi / 180.0, 1.0, 2.0);
    Matrix4d Q = Matrix4d::Zero();
    Q.diagonal() = maximum_state.array().square().inverse().matrix();
    const LqrResult lqr =
        solveCareHamiltonian(linear.A, linear.B, Q, 1.0 / 10000.0);
    RowVector4d expected;
    expected << -500.0, 4773.61252644, -605.90674606, 960.56980911;
    if ((lqr.K - expected).cwiseAbs().maxCoeff() > 1e-6) {
        std::cerr << "LQR gain regression failed\n";
        return 1;
    }
    if (lqr.care_relative_residual > 1e-10) {
        std::cerr << "CARE residual regression failed\n";
        return 1;
    }
    return 0;
}
