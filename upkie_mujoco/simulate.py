#!/usr/bin/env python3
"""Run the official Upkie model in MuJoCo with direct wheel torques."""

from __future__ import annotations

import argparse
import time

import mujoco
import numpy as np

from model import WHEEL_TORQUE_LIMIT, build_model, set_neutral_leg_targets

TORQUE_STEP = 0.05


def apply_wheel_torques(
    model: mujoco.MjModel,
    data: mujoco.MjData,
    longitudinal_torque: float,
    turn_torque: float,
) -> None:
    """Apply Upkie wheel torques in forward/turn coordinates."""
    left = np.clip(
        longitudinal_torque + turn_torque,
        -WHEEL_TORQUE_LIMIT,
        WHEEL_TORQUE_LIMIT,
    )
    right = np.clip(
        -longitudinal_torque + turn_torque,
        -WHEEL_TORQUE_LIMIT,
        WHEEL_TORQUE_LIMIT,
    )
    data.ctrl[model.actuator("left_wheel_torque").id] = left
    data.ctrl[model.actuator("right_wheel_torque").id] = right


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--headless", action="store_true", help="뷰어 없이 실행")
    parser.add_argument(
        "--duration",
        type=float,
        default=5.0,
        help="실행 시간(초). 뷰어에서 0이면 닫을 때까지 실행",
    )
    parser.add_argument(
        "--torque",
        type=float,
        default=0.0,
        help="초기 전진 토크(N·m): 왼쪽 +, 오른쪽 -",
    )
    parser.add_argument(
        "--turn-torque",
        type=float,
        default=0.0,
        help="초기 회전 토크(N·m): 양쪽 같은 부호",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model = build_model()
    data = mujoco.MjData(model)
    set_neutral_leg_targets(model, data)
    longitudinal_torque = args.torque
    turn_torque = args.turn_torque
    apply_wheel_torques(
        model, data, longitudinal_torque, turn_torque
    )
    mujoco.mj_forward(model, data)

    if args.headless:
        if args.duration <= 0.0:
            raise SystemExit("--headless에서는 --duration을 0보다 크게 지정해야 합니다.")
        while data.time < args.duration:
            mujoco.mj_step(model, data)
    else:
        from mujoco import viewer as mj_viewer

        def key_callback(keycode: int) -> None:
            nonlocal longitudinal_torque, turn_torque
            key = chr(keycode).upper() if 0 <= keycode < 256 else ""
            if key == "W":
                longitudinal_torque += TORQUE_STEP
            elif key == "S":
                longitudinal_torque -= TORQUE_STEP
            elif key == "A":
                turn_torque -= TORQUE_STEP
            elif key == "D":
                turn_torque += TORQUE_STEP
            elif keycode == 32:
                longitudinal_torque = turn_torque = 0.0
            apply_wheel_torques(
                model, data, longitudinal_torque, turn_torque
            )
            print(
                f"\rforward={longitudinal_torque:+.2f} N·m, "
                f"turn={turn_torque:+.2f} N·m",
                end="",
                flush=True,
            )

        print("W/S: 전후 토크 | A/D: 회전 토크 | Space: 토크 0")
        start_time = data.time
        with mj_viewer.launch_passive(
            model, data, key_callback=key_callback
        ) as viewer:
            viewer.cam.type = mujoco.mjtCamera.mjCAMERA_TRACKING
            viewer.cam.trackbodyid = model.body("base").id
            viewer.cam.distance = 1.5
            viewer.cam.azimuth = 135
            viewer.cam.elevation = -20
            while viewer.is_running():
                step_start = time.perf_counter()
                if args.duration > 0.0 and data.time - start_time >= args.duration:
                    break
                mujoco.mj_step(model, data)
                viewer.sync()
                remaining = model.opt.timestep - (
                    time.perf_counter() - step_start
                )
                if remaining > 0.0:
                    time.sleep(remaining)

    rotation = data.xmat[model.body("base").id].reshape(3, 3)
    pitch = np.degrees(
        np.arctan2(
            -rotation[2, 0],
            np.hypot(rotation[0, 0], rotation[1, 0]),
        )
    )
    print(
        f"\n완료: t={data.time:.3f} s, "
        f"base_xyz={data.qpos[:3]}, pitch={pitch:+.2f} deg"
    )


if __name__ == "__main__":
    main()
