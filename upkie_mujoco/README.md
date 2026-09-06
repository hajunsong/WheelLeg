# Upkie MuJoCo LQR

공식 `upkie-description` 패키지의 URDF 질량·관성·충돌 형상을 불러와
MuJoCo에서 Upkie를 시뮬레이션하고 바퀴 LQI/LQR로 균형·주행시키는
독립 실행 예제다.

## 모델 구성

- 공식 Upkie URDF의 6개 관절 사용
  - `left_hip`, `left_knee`, `left_wheel`
  - `right_hip`, `right_knee`, `right_wheel`
- 전체 질량: 약 5.34 kg
- 힙·무릎: 0 rad 중립 자세를 position PD로 유지
- 바퀴: `±1.7 N·m` 직접 토크 입력
- 로봇 루트: 6자유도 free joint
- 바퀴와 바닥: URDF cylinder 충돌 형상과 MuJoCo plane 접촉

`model.py`는 실행 시 URDF를 `MjSpec`으로 읽고 모델을 메모리에서
컴파일한다. 이렇게 해야 URDF의 fixed link가 합쳐질 때 몸체의 질량과
관성이 모두 보존된다.

## 설치

```bash
cd /home/keti/WheelLeg/upkie_mujoco
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
```

현재 폴더에는 위 가상환경이 이미 준비되어 있다.

## 실행

LQR 정지 균형:

```bash
./.venv/bin/python lqr_drive.py
```

초기 목표 속도와 방향을 지정한 주행:

```bash
./.venv/bin/python lqr_drive.py --speed 0.3 --heading 30
```

뷰어 없이 10초 검증:

```bash
./.venv/bin/python lqr_drive.py \
  --headless --duration 10 --speed 0.3 --heading 30
```

제어기 없이 원시 MuJoCo 모델만 실행:

```bash
./.venv/bin/python simulate.py --duration 0
```

원시 모델은 바퀴 균형 제어가 없으므로 중력에 의해 넘어지는 것이 정상이다.

## 키 조작

`lqr_drive.py`:

- `W/S`: 목표 전진 속도 `±0.05 m/s`
- `A/D`: 목표 방향 `±5 deg`
- `Space`: 목표 속도 0
- `H`: 현재 방향을 목표 방향으로 지정
- `R`: 초기 2도 pitch 외란으로 재시작

`simulate.py`:

- `W/S`: 전후 바퀴 토크 증감
- `A/D`: 차동 회전 토크 증감
- `Space`: 바퀴 토크 0

## 제어 구조

실행할 때 바퀴가 바닥에 닿는 직립 높이를 먼저 구한 뒤 MuJoCo의
`mjd_transitionFD`로 접촉 동역학을 수치 선형화한다.

- 전진 제어: pitch, pitch rate, 전진속도, 속도오차 적분을 사용하는 LQI
- 방향 제어: heading, yaw rate를 사용하는 LQR
- 다리 제어: 힙·무릎의 독립 position PD

Upkie의 좌우 바퀴 관절 축은 서로 반대 방향이다. 따라서 전진 토크는
`왼쪽 + / 오른쪽 -`, 회전 토크는 양쪽에 같은 부호로 적용한다.
