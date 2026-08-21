#include "lqr.hpp"
#include "model.hpp"
#include "rk4.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    using namespace wheelleg;

    try {
        // ======================= parameter =======================
        const Parameters prm = createParameters();
        constexpr double y_ref0 = 0.0;
        constexpr double v_ref = -1.0;       // cart 목표 속도 [m/s]
        constexpr double ue = 0.0;
        constexpr double force_limit = 500.0;
        const State4 xe = uprightEquilibrium(prm, y_ref0);

        const auto reference = [&xe, y_ref0, v_ref](double t) {
            State4 x_ref;
            x_ref << y_ref0 + v_ref * t, xe(1), v_ref, 0.0;
            return x_ref;
        };
        const auto autonomous_dynamics = [&prm](
                                             const State4& x,
                                             double force) {
            return nonlinearDynamics(0.0, x, prm, force);
        };

        // ======================= A, B linearization =======================
        const State4 equilibrium_residual = autonomous_dynamics(xe, ue);
        const LinearModel linear =
            numericalLinearization(autonomous_dynamics, xe, ue);
        const Eigen::Matrix4d C =
            controllabilityMatrix(linear.A, linear.B);
        Eigen::FullPivLU<Eigen::Matrix4d> controllability_lu(C);
        controllability_lu.setThreshold(1e-9);
        if (controllability_lu.rank() != 4) {
            throw std::runtime_error("linearized model is not controllable");
        }
        if (!isStabilizable(linear.A, linear.B)) {
            throw std::runtime_error("linearized model is not stabilizable");
        }

        // ======================= Q, R =======================
        constexpr double pi = 3.14159265358979323846;
        const Eigen::Vector4d maximum_state(
            0.20, 5.0 * pi / 180.0, 1.0, 2.0);
        Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
        Q.diagonal() = maximum_state.array().square().inverse().matrix();
        constexpr double force_scale = 100.0;
        constexpr double R = 1.0 / (force_scale * force_scale);

        // ======================= LQR gain =======================
        const LqrResult lqr =
            solveCareHamiltonian(linear.A, linear.B, Q, R);
        const Eigen::Matrix4d closed_loop =
            linear.A - linear.B * lqr.K;
        const Eigen::ComplexEigenSolver<Eigen::Matrix4d> pole_solver(
            closed_loop);

        std::cout << std::setprecision(10);
        std::cout << "선형화 직립 평형점 xe = "
                  << xe.transpose() << '\n';
        std::cout << "평형점 잔차 = "
                  << equilibrium_residual.transpose() << "\n\n";
        std::cout << "A =\n" << linear.A << "\n\n";
        std::cout << "B =\n" << linear.B << "\n\n";
        std::cout << "controllability rank = "
                  << controllability_lu.rank() << "/4\n";
        std::cout << "K =\n" << lqr.K << "\n\n";
        std::cout << "eig(A-BK) =\n"
                  << pole_solver.eigenvalues() << "\n";
        std::cout << "CARE relative residual = "
                  << lqr.care_relative_residual << '\n';
        std::cout << "U1 condition number = "
                  << lqr.u1_condition_number << '\n';

        const auto feedback_force =
            [&reference, &lqr, ue, force_limit](double t, const State4& x) {
            const State4 x_ref = reference(t);
            State4 error = x - x_ref;
            error(1) = angleError(x(1), x_ref(1));
            const double command = ue - (lqr.K * error)(0);
            return std::clamp(command, -force_limit, force_limit);
        };
        const std::function<State4(double, const State4&)> closed_loop_rhs =
            [&prm, &feedback_force](double t, const State4& x) {
                return nonlinearDynamics(t, x, prm, feedback_force(t, x));
            };

        // ======================= nonlinear closed-loop RK4 =======================
        constexpr double h = 0.001;
        constexpr double t_e = 10.0;
        const int n = static_cast<int>(std::lround(t_e / h));
        std::vector<double> T(n + 1);
        std::vector<State4> X(n + 1);
        std::vector<State4> Xref(n + 1);
        std::vector<double> force(n + 1);

        // 직립 상태에서 cart는 정지, pendulum 초기 각속도는 0.5 rad/s.
        State4 x = xe;
        x(3) = 0.5;
        X[0] = x;

        for (int k = 0; k <= n; ++k) {
            const double t = k * h;
            T[k] = t;
            X[k] = x;
            Xref[k] = reference(t);
            force[k] = feedback_force(t, x);
            if (k < n) {
                x = rk4Step(t, x, h, closed_loop_rhs);
            }
        }

        State4 final_error = X.back() - Xref.back();
        final_error(1) = angleError(X.back()(1), Xref.back()(1));
        double maximum_force = 0.0;
        for (const double value : force) {
            maximum_force = std::max(maximum_force, std::abs(value));
        }
        std::cout << "\n최종 추종오차 = "
                  << final_error.transpose() << '\n';
        std::cout << "최대 제어력 = " << maximum_force << " N\n";

        // ======================= CSV =======================
        const std::filesystem::path output =
            std::filesystem::path(WHEELLEG_CPP_DIR) / "lqr_response.csv";
        std::ofstream csv(output);
        csv << std::setprecision(17);
        csv << "time,y,q1,dy,dq1,y_ref,q1_ref,dy_ref,dq1_ref,force\n";
        for (int k = 0; k <= n; ++k) {
            csv << T[k];
            for (int i = 0; i < 4; ++i) {
                csv << ',' << X[k](i);
            }
            for (int i = 0; i < 4; ++i) {
                csv << ',' << Xref[k](i);
            }
            csv << ',' << force[k] << '\n';
        }
        std::cout << "LQR CSV: " << output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "오류: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
