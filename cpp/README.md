# WheelLeg C++/Eigen

MATLAB/Python 구현을 C++17과 Eigen3로 이식한 코드입니다. 모든 내부 계산은
SI 단위계(m, kg, s, N, kg·m²)를 사용합니다.

## 구성

- `main.cpp`: MATLAB `main.m` 형식의 open-loop RK4 및 RecurDyn 비교
- `dYdt.cpp`: Newton–Euler 비선형 운동방정식
- `lqr_control.cpp`: 수치 선형화, CARE, LQR 등속 추종, 비선형 RK4
- `lqr.cpp`: Eigen Hamiltonian 기반 검증형 CARE 솔버
- `model.cpp`: SI 파라미터 및 4/15차원 상태 변환
- `math_utils.cpp`: `tilde`, `ang2mat`, `mat2ep`, `ep2mat`, `step5`
- `visualize_results.py`: C++ CSV를 PNG/GIF로 변환

## 빌드

Ubuntu에서 Eigen을 시스템에 설치하려면:

```bash
sudo apt install cmake g++ libeigen3-dev
```

Eigen이 없으면 CMake가 Eigen 3.4.0을 `build/_deps`에 자동으로 내려받습니다.

```bash
cd /home/keti/WheelLeg/cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## 실행

```bash
./build/wheelleg_open_loop
./build/wheelleg_lqr
```

생성 파일:

- `open_loop.csv`: open-loop 위치·속도·가속도
- `lqr_response.csv`: LQR 실제 상태·기준 상태·제어력

LQR 목표는 `lqr_control.cpp` 상단의 다음 값으로 설정합니다.

```cpp
constexpr double y_ref0 = 0.0;
constexpr double v_ref = -1.0;
```

기준궤적은 `y_ref(t) = y_ref0 + v_ref*t`이며 진자는 직립 자세, 각속도는
0을 추종합니다.

## 시각화

Python 환경은 상위 `python/setup.sh`로 준비할 수 있습니다.

```bash
../python/.venv-linux/bin/python visualize_results.py
```

출력:

- `cpp_open_loop.png`
- `cpp_lqr_response.png`
- `cpp_lqr_animation.gif`

GIF가 필요하지 않으면 `--skip-gif`를 사용합니다.

## 검증

```bash
ctest --test-dir build --output-on-failure
```

CARE 솔버는 안정 Hamiltonian 고유공간으로 해를 계산한 뒤 다음 조건을
검증합니다.

- 안정 Hamiltonian 고유값이 정확히 4개
- 복소 허수부 및 비대칭 오차
- invariant subspace 조건수
- CARE 상대 잔차
- `P`의 반양의 정부호성
- `eig(A-BK)`의 안정성

현재 기준 이득은 다음 값과 일치해야 합니다.

```text
K = [-500, 4773.612526, -605.906746, 960.569809]
```
