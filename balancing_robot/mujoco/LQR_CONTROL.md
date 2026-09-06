# MuJoCo 밸런싱 로봇 LQR 주행 제어

## 1. 목적과 구현 범위

`lqr_drive.py`는 두 바퀴 로봇이 수직 자세를 유지하면서 목표 전진속도와
목표 heading을 추종하도록 한다. 제어기는 다음 두 채널로 분리된다.

1. **종방향 LQI**: pitch, pitch rate, 전진속도, 속도오차 적분값으로
   양쪽 바퀴의 공통 토크를 계산한다.
2. **방향 LQR**: heading 오차와 yaw rate로 좌우 바퀴의 차동 토크를 계산한다.

최종 바퀴 토크는

$$
\tau_R = u_c + u_d,\qquad
\tau_L = u_c - u_d
$$

로 구성한다. 여기서 $u_c$는 각 바퀴에 동일하게 적용되는 공통 토크,
$u_d$는 회전을 만드는 차동 토크이며 단위는 N·m이다.

## 2. 모델 및 좌표계

MuJoCo 월드 좌표계는 X=전진, Y=좌우, Z=높이이다. 원본 SolidWorks
모델은 X=전진, Y=차체 높이, Z=바퀴 축이므로 차체를 X축 기준 -90도로
회전해 Z-up 월드에 배치했다.

원본 URDF의 차체 무게중심

$$
(x,y,z)=(0.002528,\ 0,\ 0.053438)\ {\rm m}
$$

은 STL에 적용된 축 변환이 관성 좌표에는 반영되지 않은 값이었다. 이를

$$
(x,y,z)=(0.002528,\ -0.053438,\ 0)\ {\rm m}
$$

으로 수정했고, 관성 텐서에도 같은 좌표 변환을 적용했다. 제어 연습용 MJCF는
기하학적 수직 자세가 정확한 평형점이 되도록 전후 무게중심을 $x=0$으로
대칭화했다. CAD 값 2.528 mm를 그대로 사용하면 자연 평형 pitch가 약

$$
\theta_{\rm eq}
=-\tan^{-1}\left(\frac{0.002528}{0.053438}\right)
\approx -2.71^\circ
$$

가 되어, 기하학적 pitch 0도를 일정속도에서 유지하는 조건과 맞지 않는다.

바퀴 접촉은 STL 대신 반지름 0.025 m인 cylinder geom을 사용한다. 바퀴
관절은 무한 회전 hinge이며 actuator의 `gear=1`이므로 `data.ctrl` 값이
관절 토크 N·m로 직접 적용된다.

## 3. 상태 측정

차체 회전행렬을 $R_b$라 하고, 차체 로컬 축을 이용해 전방, 위쪽, 바퀴축
벡터를 다음과 같이 정의한다.

$$
\mathbf f=R_b\mathbf e_x,\qquad
\mathbf u=-R_b\mathbf e_y,\qquad
\mathbf a=R_b\mathbf e_z
$$

수평 전방 단위벡터는

$$
\mathbf f_h=
\frac{(f_x,f_y,0)^T}{\sqrt{f_x^2+f_y^2}}
$$

이고, pitch와 heading은 quaternion의 Euler angle을 직접 분해하지 않고
다음 기하식으로 계산한다.

$$
\theta=\operatorname{atan2}(\mathbf u^T\mathbf f_h,\ u_z)
$$

$$
\psi=\operatorname{atan2}(f_y,f_x)
$$

MuJoCo에서 얻은 월드 각속도 $\boldsymbol\omega$와 선속도 $\mathbf v_b$를
사용하면

$$
\dot\theta=\boldsymbol\omega^T\mathbf a,\qquad
\dot\psi=\omega_z,\qquad
v=\mathbf v_b^T\mathbf f_h
$$

이다. 이 방식은 로봇 heading이 바뀌어도 차체 기준 전진속도와 pitch를
일관되게 측정한다.

## 4. MuJoCo 수치 선형화

제어 gain은 하드코딩된 해석 모델이 아니라 실행 시 현재 MJCF에서 계산한다.
먼저 차체 회전을 수직으로 고정한 채 바퀴와 바닥의 접촉 높이를 안정화하고,
MuJoCo의 finite-difference transition 계산으로 이산 선형 모델을 얻는다.

$$
\delta\mathbf x_{k+1}
=A\,\delta\mathbf x_k+B\,\delta\mathbf u_k
$$

시뮬레이션 시간 간격은

$$
\Delta t=0.002\ {\rm s}
$$

이다. 선형화 시에는 `mjd_transitionFD`가 RK4를 지원하지 않으므로 별도로
불러온 모델의 적분기만 Euler로 바꾼다. 실제 폐루프 시뮬레이션 모델은
계속 RK4를 사용한다.

종방향 축약 상태와 공통 입력 행렬은

$$
\mathbf x_\ell=
\begin{bmatrix}
\theta & \dot\theta & v
\end{bmatrix}^T
$$

$$
B_c=B_R+B_L
$$

로 추출한다. $B_c$를 합산하는 이유는 공통 제어값 $u_c$가 양쪽 actuator에
동시에 입력되기 때문이다.

## 5. 속도·균형 LQI

목표속도 $v_{\rm ref}$에 대한 정상상태 오차를 제거하기 위해 적분 상태
$\eta$를 추가한다.

$$
e_v=v-v_{\rm ref}
$$

$$
\eta_{k+1}=\eta_k+\Delta t\,e_{v,k}
$$

증강 상태는

$$
\mathbf x_a=
\begin{bmatrix}
\theta & \dot\theta & e_v & \eta
\end{bmatrix}^T
$$

이며 선형 시스템은

$$
\mathbf x_{a,k+1}
=
\underbrace{\begin{bmatrix}
A_\ell & 0\\
\begin{matrix}0&0&\Delta t\end{matrix} & 1
\end{bmatrix}}_{A_a}
\mathbf x_{a,k}
+
\underbrace{\begin{bmatrix}B_c\\0\end{bmatrix}}_{B_a}u_{c,k}
$$

로 구성한다. 다음 이산 무한시간 비용함수를 최소화한다.

$$
J_\ell=
\sum_{k=0}^{\infty}
\left(
\mathbf x_{a,k}^TQ_\ell\mathbf x_{a,k}
+u_{c,k}^TR_\ell u_{c,k}
\right)
$$

현재 가중치는

$$
Q_\ell=\operatorname{diag}(100,\ 1,\ 5,\ 30),
\qquad R_\ell=[1]
$$

이다. 이산 Riccati 방정식의 해 $P$로 gain을 계산한다.

$$
K_\ell=
\left(R_\ell+B_a^TPB_a\right)^{-1}B_a^TPA_a
$$

$$
u_c=-K_\ell\mathbf x_a
$$

현재 모델에서 실행 시 계산되는 대표 gain은 대략

$$
K_\ell=
\begin{bmatrix}
-0.8602 & -0.07669 & -0.39568 & -0.44169
\end{bmatrix}
$$

이다. 부호는 MuJoCo 관절축과 pitch 정의에 따라 나타난 결과이며,
코드에서 다시 $u_c=-K_\ell\mathbf x_a$를 적용한다.

## 6. 방향·조향 LQR

heading 오차는 $\pm\pi$ 범위로 정규화한다.

$$
e_\psi=\operatorname{wrap}_{[-\pi,\pi)}(\psi-\psi_{\rm ref})
$$

방향 제어 상태는

$$
\mathbf x_\psi=
\begin{bmatrix}e_\psi&\dot\psi\end{bmatrix}^T
$$

이고, 차동 입력 방향은 오른쪽 $+u_d$, 왼쪽 $-u_d$로 선형화한다.

$$
\mathbf x_{\psi,k+1}
=A_\psi\mathbf x_{\psi,k}+B_\psi u_{d,k}
$$

비용함수와 가중치는

$$
J_\psi=
\sum_{k=0}^{\infty}
\left(
\mathbf x_{\psi,k}^TQ_\psi\mathbf x_{\psi,k}
+u_{d,k}^TR_\psi u_{d,k}
\right)
$$

$$
Q_\psi=\operatorname{diag}(20,\ 1),
\qquad R_\psi=[1]
$$

이다. 제어 법칙은

$$
u_d=-K_\psi\mathbf x_\psi
$$

이고 대표 gain은

$$
K_\psi=
\begin{bmatrix}0.56044&0.06509\end{bmatrix}
$$

이다.

## 7. 명령 필터와 제한

급격한 속도·방향 명령이 큰 자세 외란을 만들지 않도록 내부 기준값을 제한된
변화율로 이동한다.

$$
\left|\dot v_{\rm ref}\right|\le 0.5\ {\rm m/s^2}
$$

$$
\left|\dot\psi_{\rm ref}\right|\le 30^\circ/{\rm s}
$$

속도 적분값과 토크 제한은

$$
-1\le\eta\le1
$$

$$
|u_c|\le0.30\ {\rm N\,m},\qquad
|u_d|\le0.12\ {\rm N\,m}
$$

이다. 최종 actuator 자체에는 $[-1,1]$ N·m 제한이 추가로 적용된다.
pitch 절댓값이 60도를 넘으면 이미 넘어졌다고 판단해 바퀴 토크를 0으로 만든다.

가속하려면 접촉점에 수평 힘을 만들어야 하므로 과도구간에서 pitch가 일시적으로
0이 아닌 것은 정상이다. 목표속도에 도달한 정상상태에서는 속도 적분기가
마찰·관절 damping을 보상하고 pitch가 0도로 돌아온다.

## 8. 실행 및 조작

```bash
cd /home/keti/WheelLeg/balancing_robot/mujoco
./.venv/bin/python lqr_drive.py
```

- `W/S`: 목표속도 ±0.05 m/s
- `A/D`: 목표 heading ±5도
- `Space`: 목표속도 0
- `H`: 현재 heading 유지
- `R`: 초기 2도 pitch 외란으로 리셋

목표값을 명령행에서 지정할 수도 있다.

```bash
./.venv/bin/python lqr_drive.py --speed 0.3 --heading 45
```

뷰어 없는 회귀검증:

```bash
./.venv/bin/python lqr_drive.py \
  --headless --duration 10 --speed 0.3 --heading 45
```

검증 결과 최종 속도 0.300 m/s, heading 45.00도, pitch 약 0도로 수렴했다.
후진 -0.2 m/s와 heading -60도 조합도 같은 방식으로 추종되는 것을 확인했다.

## 9. 관련 파일

- `balancing_robot.xml`: MuJoCo 동역학, 접촉, actuator, sensor
- `lqr_drive.py`: 선형화, LQI/LQR gain 계산, 상태 측정, 실시간 제어
- `simulate.py`: 제어기 없는 수동 토크 입력 예제
- `../urdf/balancing_robot.urdf`: 수정된 원본 URDF 물성 및 연속 바퀴 관절
