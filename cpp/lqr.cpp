#include "lqr.hpp"

#include "dYdt.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

namespace wheelleg {
namespace {

double careResidual(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B,
    const Matrix4d& Q,
    double R,
    const Matrix4d& P) {
    const Matrix4d G = (B * B.transpose()) / R;
    const Matrix4d AP = A.transpose() * P;
    const Matrix4d PA = P * A;
    const Matrix4d PGP = P * G * P;
    const Matrix4d F = AP + PA - PGP + Q;
    const double denominator =
        AP.norm() + PA.norm() + PGP.norm() + Q.norm();
    return F.norm() / std::max(denominator, 1.0);
}

Matrix4d solveLyapunovCorrection(
    const Matrix4d& Ac,
    const Matrix4d& right_hand_side) {
    Eigen::Matrix<double, 16, 16> operator_matrix;
    for (int column = 0; column < 16; ++column) {
        Matrix4d basis = Matrix4d::Zero();
        const int row_index = column % 4;
        const int col_index = column / 4;
        basis(row_index, col_index) = 1.0;
        const Matrix4d mapped = Ac.transpose() * basis + basis * Ac;
        operator_matrix.col(column) =
            Eigen::Map<const Eigen::Matrix<double, 16, 1>>(mapped.data());
    }
    const Eigen::Matrix<double, 16, 1> rhs =
        Eigen::Map<const Eigen::Matrix<double, 16, 1>>(
            right_hand_side.data());
    const Eigen::Matrix<double, 16, 1> solution =
        operator_matrix.fullPivLu().solve(rhs);
    return Eigen::Map<const Matrix4d>(solution.data());
}

}  // namespace

State4 nonlinearDynamics(
    double t,
    const State4& x,
    const Parameters& prm,
    double force) {
    const State15 Yp = dYdt(t, fullState(x), prm, force);
    State4 dx;
    dx << Yp(1), Yp(7), Yp(9), Yp(14);
    return dx;
}

State4 uprightEquilibrium(const Parameters& prm, double y_ref) {
    const double q1e = std::atan2(prm.rho1p(0), prm.rho1p(1));
    State4 xe;
    xe << y_ref, q1e, 0.0, 0.0;
    return xe;
}

LinearModel numericalLinearization(
    const std::function<State4(const State4&, double)>& dynamics,
    const State4& xe,
    double ue) {
    const State4 perturbation(
        1e-6, 1e-7, 1e-6, 1e-7);
    LinearModel model;
    for (int i = 0; i < 4; ++i) {
        State4 xp = xe;
        State4 xm = xe;
        xp(i) += perturbation(i);
        xm(i) -= perturbation(i);
        model.A.col(i) =
            (dynamics(xp, ue) - dynamics(xm, ue)) /
            (2.0 * perturbation(i));
    }
    constexpr double du = 1e-4;
    model.B =
        (dynamics(xe, ue + du) - dynamics(xe, ue - du)) / (2.0 * du);
    return model;
}

Eigen::Matrix4d controllabilityMatrix(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B) {
    Eigen::Matrix4d C;
    C.col(0) = B;
    for (int i = 1; i < 4; ++i) {
        C.col(i) = A * C.col(i - 1);
    }
    return C;
}

bool isStabilizable(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B,
    double tolerance) {
    const Eigen::ComplexEigenSolver<Matrix4d> eigen_solver(A);
    if (eigen_solver.info() != Eigen::Success) {
        return false;
    }
    for (const std::complex<double>& pole : eigen_solver.eigenvalues()) {
        if (pole.real() >= -tolerance) {
            Eigen::Matrix<std::complex<double>, 4, 5> pbh;
            pbh.leftCols<4>() =
                pole * Eigen::Matrix4cd::Identity() - A.cast<std::complex<double>>();
            pbh.rightCols<1>() = B.cast<std::complex<double>>();
            Eigen::FullPivLU<Eigen::Matrix<std::complex<double>, 4, 5>> lu(pbh);
            lu.setThreshold(tolerance);
            if (lu.rank() < 4) {
                return false;
            }
        }
    }
    return true;
}

LqrResult solveCareHamiltonian(
    const Matrix4d& A,
    const Eigen::Matrix<double, 4, 1>& B,
    const Matrix4d& Q,
    double R) {
    if (!(R > 0.0) || !std::isfinite(R)) {
        throw std::runtime_error("CARE: R must be positive");
    }
    if ((Q - Q.transpose()).norm() > 1e-12 * std::max(1.0, Q.norm())) {
        throw std::runtime_error("CARE: Q must be symmetric");
    }
    const Eigen::SelfAdjointEigenSolver<Matrix4d> q_solver(Q);
    if (q_solver.info() != Eigen::Success ||
        q_solver.eigenvalues().minCoeff() < -1e-12) {
        throw std::runtime_error("CARE: Q must be positive semidefinite");
    }

    const Matrix4d G = (B * B.transpose()) / R;
    Eigen::Matrix<double, 8, 8> H;
    H.topLeftCorner<4, 4>() = A;
    H.topRightCorner<4, 4>() = -G;
    H.bottomLeftCorner<4, 4>() = -Q;
    H.bottomRightCorner<4, 4>() = -A.transpose();

    const Eigen::ComplexEigenSolver<Eigen::Matrix<double, 8, 8>> eig(H);
    if (eig.info() != Eigen::Success) {
        throw std::runtime_error("CARE: Hamiltonian eigensolver failed");
    }

    const double separation_tolerance =
        100.0 * std::sqrt(std::numeric_limits<double>::epsilon()) *
        std::max(1.0, H.norm());
    std::vector<int> stable_indices;
    for (int i = 0; i < 8; ++i) {
        const double real_part = eig.eigenvalues()(i).real();
        if (std::abs(real_part) <= separation_tolerance) {
            throw std::runtime_error(
                "CARE: Hamiltonian eigenvalue is too close to imaginary axis");
        }
        if (real_part < -separation_tolerance) {
            stable_indices.push_back(i);
        }
    }
    if (stable_indices.size() != 4) {
        throw std::runtime_error(
            "CARE: Hamiltonian does not have four stable eigenvalues");
    }

    Eigen::Matrix4cd U1;
    Eigen::Matrix4cd U2;
    for (int j = 0; j < 4; ++j) {
        U1.col(j) = eig.eigenvectors().col(stable_indices[j]).head<4>();
        U2.col(j) = eig.eigenvectors().col(stable_indices[j]).tail<4>();
    }

    const Eigen::JacobiSVD<Eigen::Matrix4cd> svd(U1);
    const auto singular = svd.singularValues();
    const double condition =
        singular(0) / std::max(singular(3), std::numeric_limits<double>::min());
    if (!std::isfinite(condition) ||
        condition * std::numeric_limits<double>::epsilon() > 1e-8) {
        throw std::runtime_error("CARE: stable invariant subspace is ill-conditioned");
    }

    // Pc = U2*U1^-1. 명시적 inverse 대신 전치 선형계를 푼다.
    const Eigen::Matrix4cd Pc =
        U1.transpose().fullPivLu().solve(U2.transpose()).transpose();
    const double imaginary_ratio =
        Pc.imag().norm() / std::max(Pc.norm(), 1.0);
    if (imaginary_ratio > 1e-8) {
        throw std::runtime_error("CARE: complex Riccati solution is not nearly real");
    }

    Matrix4d Praw = Pc.real();
    const double symmetry_ratio =
        (Praw - Praw.transpose()).norm() / std::max(Praw.norm(), 1.0);
    if (symmetry_ratio > 1e-8) {
        throw std::runtime_error("CARE: Riccati solution is not nearly symmetric");
    }
    Matrix4d P = 0.5 * (Praw + Praw.transpose());

    // Newton-Kleinman 보정: Hamiltonian 후보의 작은 수치 잔차를 제거한다.
    for (int iteration = 0; iteration < 5; ++iteration) {
        const Matrix4d F =
            A.transpose() * P + P * A - P * G * P + Q;
        if (careResidual(A, B, Q, R, P) < 1e-12) {
            break;
        }
        const Matrix4d Ac = A - G * P;
        Matrix4d correction = solveLyapunovCorrection(Ac, -F);
        correction = 0.5 * (correction + correction.transpose());
        P += correction;
    }

    const double residual = careResidual(A, B, Q, R, P);
    if (!std::isfinite(residual) || residual > 1e-9) {
        throw std::runtime_error("CARE: residual validation failed");
    }
    const Eigen::SelfAdjointEigenSolver<Matrix4d> p_solver(P);
    if (p_solver.info() != Eigen::Success ||
        p_solver.eigenvalues().minCoeff() <
            -1e-10 * std::max(1.0, p_solver.eigenvalues().maxCoeff())) {
        throw std::runtime_error("CARE: P is not positive semidefinite");
    }

    const RowVector4d K = (B.transpose() * P) / R;
    const Matrix4d closed_loop = A - B * K;
    const Eigen::ComplexEigenSolver<Matrix4d> closed_loop_eigen(closed_loop);
    if (closed_loop_eigen.info() != Eigen::Success) {
        throw std::runtime_error("CARE: closed-loop eigensolver failed");
    }
    for (const std::complex<double>& pole : closed_loop_eigen.eigenvalues()) {
        if (!(pole.real() < -1e-9)) {
            throw std::runtime_error("CARE: computed gain is not stabilizing");
        }
    }

    return {P, K, residual, condition, imaginary_ratio, symmetry_ratio};
}

double angleError(double q, double q_ref) {
    return std::atan2(std::sin(q - q_ref), std::cos(q - q_ref));
}

}  // namespace wheelleg
