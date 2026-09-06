#!/usr/bin/env python3
"""LQR balancing with velocity and heading commands."""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import mujoco
import numpy as np
from scipy.linalg import solve_discrete_are


MODEL_PATH = Path(__file__).with_name("balancing_robot.xml")
LONGITUDINAL_TORQUE_LIMIT = 0.30
TURN_TORQUE_LIMIT = 0.12
SPEED_STEP = 0.05
HEADING_STEP = math.radians(5.0)


def wrap_angle(angle: float) -> float:
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def move_toward(value: float, target: float, max_delta: float) -> float:
    return value + float(np.clip(target - value, -max_delta, max_delta))


def settle_upright(model: mujoco.MjModel) -> tuple[mujoco.MjData, float]:
    """차체 회전을 고정한 채 바퀴-바닥 접촉 높이만 안정화한다."""
    data = mujoco.MjData(model)
    reference_qpos = data.qpos.copy()
    for _ in range(1000):
        mujoco.mj_step(model, data)
        height, vertical_speed = data.qpos[2], data.qvel[2]
        data.qpos[:] = reference_qpos
        data.qpos[2] = height
        data.qvel[:] = 0.0
        data.qvel[2] = vertical_speed
    data.qvel[:] = 0.0
    mujoco.mj_forward(model, data)
    return data, float(data.qpos[2])


def calculate_lqr_gains() -> tuple[np.ndarray, np.ndarray, float]:
    """MuJoCo 접촉 동역학을 수치 선형화해 LQI/LQR gain을 계산한다."""
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
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

    # MuJoCo tangent state: pitch=q[5], pitch rate=v[5], world forward speed=v[0].
    long_indices = [5, model.nv + 5, model.nv]
    long_a = transition_a[np.ix_(long_indices, long_indices)]
    long_b = transition_b[long_indices].sum(axis=1, keepdims=True)

    # 속도 정상상태 오차 제거를 위해 적분 상태를 추가한 LQI이다.
    dt = model.opt.timestep
    augmented_a = np.block(
        [
            [long_a, np.zeros((3, 1))],
            [np.array([[0.0, 0.0, dt]]), np.ones((1, 1))],
        ]
    )
    augmented_b = np.vstack([long_b, [[0.0]]])
    long_q = np.diag([100.0, 1.0, 5.0, 30.0])
    long_r = np.array([[1.0]])
    long_p = solve_discrete_are(augmented_a, augmented_b, long_q, long_r)
    long_k = np.linalg.solve(
        long_r + augmented_b.T @ long_p @ augmented_b,
        augmented_b.T @ long_p @ augmented_a,
    ).ravel()

    # local -Y 회전이 world yaw이고, 오른쪽(+)/왼쪽(-) 토크가 양의 yaw이다.
    yaw_indices = [4, model.nv + 4]
    yaw_a = transition_a[np.ix_(yaw_indices, yaw_indices)]
    yaw_b = -(
        transition_b[yaw_indices, 0] - transition_b[yaw_indices, 1]
    )[:, None]
    yaw_q = np.diag([20.0, 1.0])
    yaw_r = np.array([[1.0]])
    yaw_p = solve_discrete_are(yaw_a, yaw_b, yaw_q, yaw_r)
    yaw_k = np.linalg.solve(
        yaw_r + yaw_b.T @ yaw_p @ yaw_b,
        yaw_b.T @ yaw_p @ yaw_a,
    ).ravel()
    return long_k, yaw_k, upright_height


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
        self.base_id = model.body("base_link").id
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
        up = -rotation[:, 1]
        lateral = rotation[:, 2]
        horizontal_forward = forward.copy()
        horizontal_forward[2] = 0.0
        horizontal_forward /= np.linalg.norm(horizontal_forward)

        pitch = math.atan2(
            float(up @ horizontal_forward),
            float(up[2]),
        )
        heading = math.atan2(float(forward[1]), float(forward[0]))

        spatial_velocity = np.empty(6)
        mujoco.mj_objectVelocity(
            self.model,
            data,
            mujoco.mjtObj.mjOBJ_BODY,
            self.base_id,
            spatial_velocity,
            0,
        )
        angular_velocity = spatial_velocity[:3]
        linear_velocity = spatial_velocity[3:]
        pitch_rate = float(angular_velocity @ lateral)
        yaw_rate = float(angular_velocity[2])
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

        # 급격한 명령으로 넘어지지 않도록 속도와 방향 기준값을 제한 속도로 이동한다.
        self.speed_reference = move_toward(
            self.speed_reference, self.target_speed, 0.5 * dt
        )
        heading_delta = wrap_angle(self.target_heading - self.heading_reference)
        self.heading_reference += float(
            np.clip(heading_delta, -math.radians(30.0) * dt, math.radians(30.0) * dt)
        )

        speed_error = speed - self.speed_reference
        self.speed_integral = float(
            np.clip(self.speed_integral + dt * speed_error, -1.0, 1.0)
        )
        long_state = np.array(
            [pitch, pitch_rate, speed_error, self.speed_integral]
        )
        common_torque = float(
            np.clip(
                -self.longitudinal_gain @ long_state,
                -LONGITUDINAL_TORQUE_LIMIT,
                LONGITUDINAL_TORQUE_LIMIT,
            )
        )

        yaw_state = np.array(
            [wrap_angle(heading - self.heading_reference), yaw_rate]
        )
        turn_torque = float(
            np.clip(
                -self.yaw_gain @ yaw_state,
                -TURN_TORQUE_LIMIT,
                TURN_TORQUE_LIMIT,
            )
        )
        right_torque = common_torque + turn_torque
        left_torque = common_torque - turn_torque
        return right_torque, left_torque, state


def reset_robot(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    upright_height: float,
    initial_tilt: float,
) -> None:
    mujoco.mj_resetData(model, data)
    data.qpos[2] = upright_height
    angle = math.radians(initial_tilt)
    tilt_quaternion = np.array(
        [math.cos(angle / 2.0), 0.0, math.sin(angle / 2.0), 0.0]
    )
    result = np.empty(4)
    mujoco.mju_mulQuat(result, tilt_quaternion, data.qpos[3:7])
    data.qpos[3:7] = result
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

    print("수직 평형점에서 MuJoCo 동역학을 선형화하고 있습니다...")
    longitudinal_gain, yaw_gain, upright_height = calculate_lqr_gains()
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
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
    right_id = model.actuator("right_wheel_torque").id
    left_id = model.actuator("left_wheel_torque").id

    def step_control() -> tuple[float, float, float, float, float]:
        right_torque, left_torque, state = controller.control(data)
        if abs(state[0]) > math.radians(60.0):
            right_torque = left_torque = 0.0
        data.ctrl[right_id] = right_torque
        data.ctrl[left_id] = left_torque
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
            viewer.cam.distance = 0.6
            viewer.cam.azimuth = 135
            viewer.cam.elevation = -20
            while viewer.is_running():
                step_start = time.perf_counter()
                if args.duration > 0.0 and data.time - start_time >= args.duration:
                    break
                state = step_control()
                if data.time >= next_status_time:
                    print(
                        f"\rspeed={state[2]:+.3f}/{controller.target_speed:+.3f} m/s, "
                        f"heading={math.degrees(state[3]):+.1f}/"
                        f"{math.degrees(controller.target_heading):+.1f} deg, "
                        f"pitch={math.degrees(state[0]):+.2f} deg",
                        end="",
                        flush=True,
                    )
                    next_status_time = data.time + 0.2
                viewer.sync()
                remaining = model.opt.timestep - (time.perf_counter() - step_start)
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
