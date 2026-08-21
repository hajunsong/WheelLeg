#pragma once

#include "model.hpp"

#include <functional>

namespace wheelleg {

inline State15 rk4Step(
    double t,
    const State15& Y,
    double h,
    const std::function<State15(double, const State15&)>& rhs) {
    const State15 k1 = rhs(t, Y);
    const State15 k2 = rhs(t + h / 2.0, Y + h * k1 / 2.0);
    const State15 k3 = rhs(t + h / 2.0, Y + h * k2 / 2.0);
    const State15 k4 = rhs(t + h, Y + h * k3);
    State15 next = Y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    next.segment<4>(3).normalize();
    return next;
}

inline State4 rk4Step(
    double t,
    const State4& x,
    double h,
    const std::function<State4(double, const State4&)>& rhs) {
    const State4 k1 = rhs(t, x);
    const State4 k2 = rhs(t + h / 2.0, x + h * k1 / 2.0);
    const State4 k3 = rhs(t + h / 2.0, x + h * k2 / 2.0);
    const State4 k4 = rhs(t + h, x + h * k3);
    return x + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

}  // namespace wheelleg
