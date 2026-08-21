#include "dYdt.hpp"
#include "model.hpp"
#include "rk4.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Sample6 = std::array<double, 6>;

struct ReferenceData {
    std::vector<double> time;
    std::vector<Sample6> value;
};

ReferenceData readRecurDynCsv(const std::filesystem::path& path) {
    ReferenceData ref;
    std::ifstream file(path);
    if (!file) {
        return ref;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream stream(line);
        std::array<double, 8> row{};
        std::string token;
        for (double& value : row) {
            if (!std::getline(stream, token, ',')) {
                throw std::runtime_error("invalid RecurDyn CSV row");
            }
            value = std::stod(token);
        }
        ref.time.push_back(row[1]);
        ref.value.push_back({
            1e-3 * row[2], 1e-3 * row[3], 1e-3 * row[4],
            row[5], row[6], row[7],
        });
    }
    return ref;
}

double interpolate(
    const std::vector<double>& time,
    const std::vector<Sample6>& value,
    double query,
    int channel) {
    if (query <= time.front()) {
        return value.front()[channel];
    }
    if (query >= time.back()) {
        return value.back()[channel];
    }
    const auto upper = std::upper_bound(time.begin(), time.end(), query);
    const std::size_t hi = static_cast<std::size_t>(upper - time.begin());
    const std::size_t lo = hi - 1;
    const double alpha = (query - time[lo]) / (time[hi] - time[lo]);
    return (1.0 - alpha) * value[lo][channel] + alpha * value[hi][channel];
}

}  // namespace

int main() {
    using namespace wheelleg;

    try {
        // ======================= parameter =======================
        // 단위계: SI (m, kg, s, N, kg*m^2)
        const Parameters prm = createParameters();
        const double h = 0.001;
        const double t_e = 2.0;

        // ======================= initial condition =======================
        State4 x0;
        x0 << 0.0, 0.0, 0.0, 3.0;
        State15 Y = fullState(x0);

        // ======================= RK4 =======================
        const int n = static_cast<int>(std::lround(t_e / h));
        std::vector<double> T(n + 1);
        std::vector<State15> YY(n + 1);
        std::vector<State15> AA(n + 1);
        YY[0] = Y;

        const auto rhs = [&prm](double t, const State15& state) {
            return dYdt(t, state, prm);
        };
        for (int k = 0; k < n; ++k) {
            T[k] = k * h;
            AA[k] = rhs(T[k], Y);
            Y = rk4Step(T[k], Y, h, rhs);
            YY[k + 1] = Y;
        }
        T[n] = t_e;
        AA[n] = rhs(T[n], YY[n]);

        // ======================= post processing =======================
        std::vector<Sample6> mine(n + 1);
        double max_quaternion_drift = 0.0;
        for (int k = 0; k <= n; ++k) {
            mine[k] = {
                YY[k](1), YY[k](9), AA[k](9),
                YY[k](7), YY[k](14), AA[k](14),
            };
            max_quaternion_drift = std::max(
                max_quaternion_drift,
                std::abs(YY[k].segment<4>(3).norm() - 1.0));
        }

        const std::filesystem::path output =
            std::filesystem::path(WHEELLEG_CPP_DIR) / "open_loop.csv";
        std::ofstream csv(output);
        csv << std::setprecision(17);
        csv << "time,cart_y,cart_vy,cart_ay,q1,dq1,ddq1\n";
        for (int k = 0; k <= n; ++k) {
            csv << T[k];
            for (double value : mine[k]) {
                csv << ',' << value;
            }
            csv << '\n';
        }

        std::cout << std::setprecision(8);
        std::cout << "t = " << T.front() << " ~ " << T.back()
                  << " s, " << n << " steps (h = " << h << ")\n";
        std::cout << "quaternion norm drift : "
                  << max_quaternion_drift << '\n';
        std::cout << "open-loop CSV: " << output << '\n';

        // ======================= RecurDyn 비교 =======================
        const std::filesystem::path ref_path =
            std::filesystem::path(WHEELLEG_SOURCE_DIR) /
            "recurdyn/01_cart_pole/rec_data.csv";
        const ReferenceData ref = readRecurDynCsv(ref_path);
        if (!ref.time.empty()) {
            const std::array<std::string, 6> labels = {
                "cart_y [m]", "cart_vy [m/s]", "cart_ay [m/s^2]",
                "q1 [rad]", "dq1 [rad/s]", "ddq1 [rad/s^2]",
            };
            std::cout << "\nRecurDyn: " << ref.value.size() << " points\n";
            std::cout << std::left << std::setw(22) << "channel"
                      << std::right << std::setw(16) << "max|err|"
                      << std::setw(16) << "RMS"
                      << std::setw(16) << "rel.RMS" << '\n';
            for (int channel = 0; channel < 6; ++channel) {
                double max_error = 0.0;
                double sum_square = 0.0;
                double scale = 0.0;
                for (std::size_t i = 0; i < ref.time.size(); ++i) {
                    const double predicted =
                        interpolate(T, mine, ref.time[i], channel);
                    const double error = predicted - ref.value[i][channel];
                    max_error = std::max(max_error, std::abs(error));
                    sum_square += error * error;
                    scale = std::max(scale, std::abs(ref.value[i][channel]));
                }
                const double rms =
                    std::sqrt(sum_square / static_cast<double>(ref.time.size()));
                const double relative = scale > 0.0 ? rms / scale : 0.0;
                std::cout << std::left << std::setw(22) << labels[channel]
                          << std::right << std::scientific
                          << std::setw(16) << max_error
                          << std::setw(16) << rms
                          << std::setw(16) << relative << '\n';
            }
        } else {
            std::cout << "\n[비교 생략] " << ref_path << " 없음\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "오류: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
