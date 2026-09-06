#!/usr/bin/env python3
"""Balance and drive the official Upkie model with LQI/LQR control."""

from __future__ import annotations

import argparse
import math
import time

import mujoco
import numpy as np
from scipy.linalg import solve_discrete_are

from model import WHEEL_TORQUE_LIMIT, build_model, set_neutral_leg_targets

LONGITUDINAL_TORQUE_LIMIT = 1.5
TURN_TORQUE_LIMIT = 0.1
SPEED_STEP = 0.05
HEADING_STEP = math.radians(5.0)


def wrap_angle(angle: float) -> float:
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def move_toward(value: float, target: float, max_delta: float) -> float:
    return value + float(np.clip(target - value, -max_delta, max_delta))


def settle_upright(model: mujoco.MjModel) -> tuple[mujoco.MjData, float]:
    """Find wheel contact height while holding the nominal configuration."""
    data = mujoco.MjData(model)
    set_neutral_leg_targets(model, data)
    reference_qpos = data.qpos.copy()
    for _ in range(1500):
        mujoco.mj_step(model, data)
        height = data.qpos[2]
        vertical_speed = data.qvel[2]
        data.qpos[:] = reference_qpos
        data.qpos[2] = height
        data.qvel[:] = 0.0
        data.qvel[2] = vertical_speed
        set_neutral_leg_targets(model, data)
    data.qvel[:] = 0.0
    mujoco.mj_forward(model, data)
    return data, float(data.qpos[2])


def _discrete_lqr_gain(
    transition_a: np.ndarray,
    transition_b: np.ndarray,
    state_cost: np.ndarray,
    input_cost: np.ndarray,
) -> np.ndarray:
    value = solve_discrete_are(
        transition_a, transition_b, state_cost, input_cost
    )
    return np.linalg.solve(
        input_cost + transition_b.T @ value @ transition_b,
        transition_b.T @ value @ transition_a,
    ).ravel()


def calculate_lqr_gains() -> tuple[np.ndarray, np.ndarray, float]:
    """Linearize wheel-contact dynamics and calculate drive controller gains."""
    model = build_model()
    data, upright_height = settle_upright(model)
    model.opt.integrator = mujoco.mjtIntegrator.mjINT_EULER

    state_size = 2 * model.nv + model.na
    transition_a = np.empty((state_size, state_size))
    transition_b = np.empty((state_size, model.nu))
    mujoco.mjd_transitionFD(
        model,
        data,
        1e-6,
        1,
        transition_a,
        transition_b,
        None,
        None,
    )

    free_dof = model.joint("base_free_joint").dofadr[0]
    pitch_dof = free_dof + 4
    yaw_dof = free_dof + 5
    forward_dof = free_dof
    left_id = model.actuator("left_wheel_torque").id
    right_id = model.actuator("right_wheel_torque").id

    # Upkie is left-wheeled: forward torque is +left and -right.
    long_indices = [pitch_dof, model.nv + pitch_dof, model.nv + forward_dof]
    long_a = transition_a[np.ix_(long_indices, long_indices)]
    long_b = (
        transition_b[long_indices, left_id]
        - transition_b[long_indices, right_id]
    )[:, None]

    dt = model.opt.timestep
    augmented_a = np.block(
        [
            [long_a, np.zeros((3, 1))],
            [np.array([[0.0, 0.0, dt]]), np.ones((1, 1))],
        ]
    )
    augmented_b = np.vstack([long_b, [[0.0]]])
    long_gain = _discrete_lqr_gain(
        augmented_a,
        augmented_b,
        np.diag([200.0, 5.0, 10.0, 30.0]),
        np.array([[1.0]]),
    )

    # Equal rotor torques create the opposite base-yaw reaction through the
    # wheel contacts, hence the minus sign in this reduced input coordinate.
    yaw_indices = [yaw_dof, model.nv + yaw_dof]
    yaw_a = transition_a[np.ix_(yaw_indices, yaw_indices)]
    yaw_b = -(
        transition_b[yaw_indices, left_id]
        + transition_b[yaw_indices, right_id]
    )[:, None]
    yaw_gain = _discrete_lqr_gain(
        yaw_a,
        yaw_b,
        np.diag([20.0, 1.0]),
        np.array([[100.0]]),
    )
    return long_gain, yaw_gain, upright_height


class LQRDriveController:
    def __init__(
        self,
        model: mujoco.MjModel,
        longitudinal_gain: np.ndarray,
        yaw_gain: np.ndarray,
        target_speed: float,
        target_heading: float,
    ) -> None:
        self.model = model
        self.base_id = model.body("base").id
        self.left_wheel_id = model.actuator("left_wheel_torque").id
        self.right_wheel_id = model.actuator("right_wheel_torque").id
        self.free_dof = model.joint("base_free_joint").dofadr[0]
        self.longitudinal_gain = longitudinal_gain
        self.yaw_gain = yaw_gain
        self.target_speed = target_speed
        self.target_heading = target_heading
        self.speed_reference = 0.0
        self.heading_reference = 0.0
        self.speed_integral = 0.0

    def measure(
        self, data: mujoco.MjData
    ) -> tuple[float, float, float, float, float]:
        rotation = data.xmat[self.base_id].reshape(3, 3)
        forward = rotation[:, 0]
        horizontal_forward = forward.copy()
        horizontal_forward[2] = 0.0
        norm = np.linalg.norm(horizontal_forward)
        if norm > 1e-9:
            horizontal_forward /= norm

        pitch = math.atan2(
            -float(forward[2]),
            math.hypot(float(forward[0]), float(forward[1])),
        )
        heading = math.atan2(float(forward[1]), float(forward[0]))

        # Free-joint velocity ordering is translation XYZ then rotation XYZ.
        # Reading this state directly also matches the coordinates used by
        # mjd_transitionFD during controller linearization.
        linear_velocity = data.qvel[self.free_dof : self.free_dof + 3]
        pitch_rate = float(data.qvel[self.free_dof + 4])
        yaw_rate = float(data.qvel[self.free_dof + 5])
        forward_speed = float(linear_velocity @ horizontal_forward)
        return pitch, pitch_rate, forward_speed, heading, yaw_rate

    def reset(self, data: mujoco.MjData) -> None:
        _, _, _, heading, _ = self.measure(data)
        self.speed_reference = 0.0
        self.heading_reference = heading
        self.speed_integral = 0.0

    def control(
        self, data: mujoco.MjData
    ) -> tuple[float, float, tuple[float, float, float, float, float]]:
        dt = self.model.opt.timestep
        state = self.measure(data)
        pitch, pitch_rate, speed, heading, yaw_rate = state

        self.speed_reference = move_toward(
            self.speed_reference, self.target_speed, 0.5 * dt
        )
        heading_delta = wrap_angle(self.target_heading - self.heading_reference)
        self.heading_reference += float(
            np.clip(
                heading_delta,
                -math.radians(30.0) * dt,
                math.radians(30.0) * dt,
            )
        )

        speed_error = speed - self.speed_reference
        self.speed_integral = float(
            np.clip(self.speed_integral + dt * speed_error, -1.0, 1.0)
        )
        common_torque = float(
            np.clip(
                -self.longitudinal_gain
                @ np.array(
                    [pitch, pitch_rate, speed_error, self.speed_integral]
                ),
                -LONGITUDINAL_TORQUE_LIMIT,
                LONGITUDINAL_TORQUE_LIMIT,
            )
        )

        turn_torque = float(
            np.clip(
                -self.yaw_gain
                @ np.array(
                    [wrap_angle(heading - self.heading_reference), yaw_rate]
                ),
                -TURN_TORQUE_LIMIT,
                TURN_TORQUE_LIMIT,
            )
        )
        left_torque = float(
            np.clip(
                common_torque + turn_torque,
                -WHEEL_TORQUE_LIMIT,
                WHEEL_TORQUE_LIMIT,
            )
        )
        right_torque = float(
            np.clip(
                -common_torque + turn_torque,
                -WHEEL_TORQUE_LIMIT,
                WHEEL_TORQUE_LIMIT,
            )
        )
        return left_torque, right_torque, state


def reset_robot(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    upright_height: float,
    initial_tilt: float,
) -> None:
    mujoco.mj_resetData(model, data)
    data.qpos[2] = upright_height
    angle = math.radians(initial_tilt)
    data.qpos[3:7] = [
        math.cos(angle / 2.0),
        0.0,
        math.sin(angle / 2.0),
        0.0,
    ]
    set_neutral_leg_targets(model, data)
    mujoco.mj_forward(model, data)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--speed", type=float, default=0.0, help="목표 속도(m/s)")
    parser.add_argument(
        "--heading", type=float, default=0.0, help="목표 방향(deg, 반시계 +)"
    )
    parser.add_argument(
        "--initial-tilt", type=float, default=2.0, help="초기 pitch 외란(deg)"
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.0,
        help="실행 시간(초). 0이면 뷰어를 닫을 때까지 실행",
    )
    parser.add_argument("--headless", action="store_true", help="뷰어 없이 실행")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.headless and args.duration <= 0.0:
        raise SystemExit("--headless에서는 --duration을 0보다 크게 지정해야 합니다.")

    print("Upkie 직립점에서 MuJoCo 동역학을 선형화하고 있습니다...")
    longitudinal_gain, yaw_gain, upright_height = calculate_lqr_gains()
    print(f"longitudinal gain: {longitudinal_gain}")
    print(f"yaw gain: {yaw_gain}")

    model = build_model()
    data = mujoco.MjData(model)
    reset_robot(model, data, upright_height, args.initial_tilt)
    controller = LQRDriveController(
        model,
        longitudinal_gain,
        yaw_gain,
        args.speed,
        math.radians(args.heading),
    )
    controller.reset(data)

    def step_control() -> tuple[float, float, float, float, float]:
        left_torque, right_torque, state = controller.control(data)
        if abs(state[0]) > math.radians(60.0):
            left_torque = right_torque = 0.0
        data.ctrl[controller.left_wheel_id] = left_torque
        data.ctrl[controller.right_wheel_id] = right_torque
        mujoco.mj_step(model, data)
        return state

    if args.headless:
        while data.time < args.duration:
            state = step_control()
    else:
        from mujoco import viewer as mj_viewer

        def key_callback(keycode: int) -> None:
            key = chr(keycode).upper() if 0 <= keycode < 256 else ""
            if key == "W":
                controller.target_speed += SPEED_STEP
            elif key == "S":
                controller.target_speed -= SPEED_STEP
            elif key == "A":
                controller.target_heading += HEADING_STEP
            elif key == "D":
                controller.target_heading -= HEADING_STEP
            elif keycode == 32:
                controller.target_speed = 0.0
            elif key == "H":
                controller.target_heading = controller.measure(data)[3]
            elif key == "R":
                reset_robot(model, data, upright_height, args.initial_tilt)
                controller.reset(data)

        print(
            "W/S: 목표속도 ±0.05 m/s | A/D: 목표방향 ±5 deg | "
            "Space: 정지 | H: 현재 방향 유지 | R: 리셋"
        )
        start_time = data.time
        next_status_time = 0.0
        with mj_viewer.launch_passive(
            model, data, key_callback=key_callback
        ) as viewer:
            viewer.cam.type = mujoco.mjtCamera.mjCAMERA_TRACKING
            viewer.cam.trackbodyid = controller.base_id
            viewer.cam.distance = 1.5
            viewer.cam.azimuth = 135
            viewer.cam.elevation = -20
            while viewer.is_running():
                step_start = time.perf_counter()
                if args.duration > 0.0 and data.time - start_time >= args.duration:
                    break
                state = step_control()
                if data.time >= next_status_time:
                    print(
                        f"\rspeed={state[2]:+.3f}/"
                        f"{controller.target_speed:+.3f} m/s, "
                        f"heading={math.degrees(state[3]):+.1f}/"
                        f"{math.degrees(controller.target_heading):+.1f} deg, "
                        f"pitch={math.degrees(state[0]):+.2f} deg",
                        end="",
                        flush=True,
                    )
                    next_status_time = data.time + 0.2
                viewer.sync()
                remaining = model.opt.timestep - (
                    time.perf_counter() - step_start
                )
                if remaining > 0.0:
                    time.sleep(remaining)

    final_state = controller.measure(data)
    print(
        f"\n완료: speed={final_state[2]:+.3f} m/s, "
        f"heading={math.degrees(final_state[3]):+.2f} deg, "
        f"pitch={math.degrees(final_state[0]):+.3f} deg"
    )


if __name__ == "__main__":
    main()
