# MuJoCo balancing robot

SolidWorks에서 생성된 `../urdf/balancing_robot.urdf`의 질량, 관성 및 STL 메시를
사용한 MuJoCo 모델이다. 원본 모델 좌표계를 회전해 MuJoCo 월드는
X=전진, Y=좌우, Z=높이(Z-up)가 되도록 구성했다.

## 설치 및 실행

```bash
cd balancing_robot/mujoco
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt

# 뷰어 실행
./.venv/bin/python simulate.py --duration 0

# 뷰어 없이 양쪽 바퀴에 0.02 N·m를 2초 동안 적용
./.venv/bin/python simulate.py --headless --duration 2 --torque 0.02

# LQR balancing 제어(창을 닫을 때까지)
./.venv/bin/python lqr_drive.py

# 0.3 m/s, 반시계 방향 45도 목표
./.venv/bin/python lqr_drive.py --speed 0.3 --heading 45

# 동일 명령을 뷰어 없이 10초 검증
./.venv/bin/python lqr_drive.py --headless --duration 10 --speed 0.3 --heading 45
```

뷰어에서는 `W/S`로 양쪽 토크를 증감하고, `A/D`로 좌우 차동 토크를 주며,
`Space`로 토크를 0으로 만든다. 토크는 한 번 누를 때 0.02 N·m씩 바뀐다.

## LQR 주행 제어

`lqr_drive.py`는 실행할 때 MuJoCo의 수직 평형점 접촉 동역학을 수치 선형화한다.
공통 바퀴 토크에는 pitch·pitch rate·전진속도·속도 적분 상태를 사용하는 LQI를,
차동 토크에는 heading·yaw rate LQR을 적용한다.

- `W/S`: 목표 속도 `±0.05 m/s`
- `A/D`: 목표 방향 `±5 deg`
- `Space`: 목표 속도 0
- `H`: 현재 방향을 새 목표 방향으로 지정
- `R`: 초기 2도 pitch 외란 상태로 리셋

급격한 가속 중에는 바퀴의 가속력을 만들기 위한 pitch 편차가 일시적으로 발생한다.
목표 속도에 도달하면 pitch는 0도로 복귀한다.

## 제어 입력

`data.ctrl`에 입력한 값이 각 바퀴의 관절 토크(N·m)로 직접 전달된다.

```python
right_id = model.actuator("right_wheel_torque").id
left_id = model.actuator("left_wheel_torque").id
data.ctrl[right_id] = right_torque_nm
data.ctrl[left_id] = left_torque_nm
mujoco.mj_step(model, data)
```

두 actuator의 입력 범위는 `balancing_robot.xml`에서 `[-1, 1] N·m`로 제한했다.
실물 모터 사양에 맞춰 `ctrlrange`와 `TORQUE_LIMIT`을 함께 조정하면 된다.

## 제어에 사용할 수 있는 센서

- `base_orientation`: 차체 quaternion `(w, x, y, z)`
- `base_angular_velocity`: 차체 각속도 `(rad/s)`
- `base_acceleration`: IMU 위치의 선가속도
- `right_wheel_position`, `left_wheel_position`: 바퀴 각도 `(rad)`
- `right_wheel_velocity`, `left_wheel_velocity`: 바퀴 각속도 `(rad/s)`
- `right_wheel_applied_torque`, `left_wheel_applied_torque`: 적용 토크 `(N·m)`

이름으로 센서 ID를 얻거나 `data.sensor("센서명").data`로 값을 읽을 수 있다.

저장소의 기존 `python/lqr_control.py`는 이 로봇이 아니라 질량과 입력 방식이 다른
cart-pendulum 모델(`Fy` 선형 힘 입력)을 대상으로 한다. 따라서 바퀴 토크 제어기에
그 gain을 그대로 사용하면 안 되며, 이 MuJoCo 모델을 기준으로 다시 선형화해야 한다.

## 모델 구성 메모

- 바퀴는 무한 회전하는 hinge joint이다.
- STL은 시각화에 사용하고, 접촉 계산에는 box/cylinder 단순 형상을 사용한다.
- 로봇 루트는 free joint이므로 3차원 이동과 넘어짐이 가능하다.
- 원본 URDF의 바퀴 관절도 `continuous`로 수정했다.
- 원본 URDF의 차체 무게중심은 STL 축 변환이 누락돼 메시 밖에 있었으므로
  `(x, y, z)=(0.002528, -0.053438, 0)`으로 바로잡았다.
- 제어 연습용 MJCF는 수직 자세가 정확한 평형점이 되도록 전후 무게중심 X를
  0으로 대칭화했다. CAD 값 2.528 mm를 그대로 쓰면 자연 평형 tilt는 약 2.7도이다.
