#!/usr/bin/env python3
"""MuJoCo balancing robot torque-input example."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import mujoco


MODEL_PATH = Path(__file__).with_name("balancing_robot.xml")
TORQUE_LIMIT = 1.0
TORQUE_STEP = 0.02


def actuator_id(model: mujoco.MjModel, name: str) -> int:
    index = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_ACTUATOR, name)
    if index < 0:
        raise ValueError(f"actuator를 찾을 수 없습니다: {name}")
    return index


def set_wheel_torques(
    data: mujoco.MjData,
    right_id: int,
    left_id: int,
    right_torque: float,
    left_torque: float,
) -> None:
    """바퀴 토크를 N·m 단위로 적용한다."""
    data.ctrl[right_id] = max(-TORQUE_LIMIT, min(TORQUE_LIMIT, right_torque))
    data.ctrl[left_id] = max(-TORQUE_LIMIT, min(TORQUE_LIMIT, left_torque))


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
        help="양쪽 바퀴에 적용할 초기 토크(N·m)",
    )
    parser.add_argument("--right-torque", type=float, help="오른쪽 초기 토크(N·m)")
    parser.add_argument("--left-torque", type=float, help="왼쪽 초기 토크(N·m)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
    data = mujoco.MjData(model)

    right_id = actuator_id(model, "right_wheel_torque")
    left_id = actuator_id(model, "left_wheel_torque")
    right_torque = args.torque if args.right_torque is None else args.right_torque
    left_torque = args.torque if args.left_torque is None else args.left_torque
    set_wheel_torques(data, right_id, left_id, right_torque, left_torque)
    mujoco.mj_forward(model, data)

    if args.headless:
        end_time = data.time + args.duration
        while data.time < end_time:
            mujoco.mj_step(model, data)
    else:
        from mujoco import viewer as mj_viewer

        def key_callback(keycode: int) -> None:
            nonlocal right_torque, left_torque
            key = chr(keycode).upper() if 0 <= keycode < 256 else ""
            if key == "W":
                right_torque += TORQUE_STEP
                left_torque += TORQUE_STEP
            elif key == "S":
                right_torque -= TORQUE_STEP
                left_torque -= TORQUE_STEP
            elif key == "A":
                right_torque += TORQUE_STEP
                left_torque -= TORQUE_STEP
            elif key == "D":
                right_torque -= TORQUE_STEP
                left_torque += TORQUE_STEP
            elif keycode == 32:
                right_torque = left_torque = 0.0
            set_wheel_torques(
                data, right_id, left_id, right_torque, left_torque
            )
            print(
                f"\rright={data.ctrl[right_id]:+.3f} N·m, "
                f"left={data.ctrl[left_id]:+.3f} N·m",
                end="",
                flush=True,
            )

        print("키 조작: W/S=전후 토크, A/D=회전 토크, Space=토크 0")
        start_time = data.time
        with mj_viewer.launch_passive(
            model, data, key_callback=key_callback
        ) as viewer:
            viewer.cam.fixedcamid = mujoco.mj_name2id(
                model, mujoco.mjtObj.mjOBJ_CAMERA, "side"
            )
            viewer.cam.type = mujoco.mjtCamera.mjCAMERA_FIXED
            while viewer.is_running():
                step_start = time.perf_counter()
                if args.duration > 0 and data.time - start_time >= args.duration:
                    break
                mujoco.mj_step(model, data)
                viewer.sync()
                remaining = model.opt.timestep - (time.perf_counter() - step_start)
                if remaining > 0:
                    time.sleep(remaining)

    print(
        f"\n완료: t={data.time:.3f} s, "
        f"base_xyz={data.qpos[:3]}, "
        f"ctrl=[{data.ctrl[right_id]:+.3f}, {data.ctrl[left_id]:+.3f}] N·m"
    )


if __name__ == "__main__":
    main()
