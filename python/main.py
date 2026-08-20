r"""main.m 1:1 변환.

실행:  .venv\Scripts\python.exe main.py
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from ang2mat import ang2mat
from dYdt import dYdt
from mat2ep import mat2ep
from util import col, rms

# ======================= parameter =======================
# 단위계 : MMKS (mm, kg, s).  힘의 내부 단위는 kg*mm/s^2 (= mN) 이므로
#          N 로 주어진 값은 1000 을 곱해서 넣는다.
prm = SimpleNamespace()

# ---- base body (기준계 = base.Ai) ----
prm.rho0p = col(0, 0, 0)
prm.C00 = ang2mat(0, np.pi/2, 0)

prm.m0 = 111.015764646288
Ixx = 9343826.85772924;  Ixy = 0.0
Iyy = 277539.41161572;   Iyz = 0.0
Izz = 9436339.99493448;  Izx = 0.0
prm.J0p = np.array([[Ixx, Ixy, Izx],
                    [Ixy, Iyy, Iyz],
                    [Izx, Iyz, Izz]])

# ---- link body (기준계 = body.Ai) ----
prm.s01p = col(0, 0, 0)
prm.C01 = ang2mat(np.pi/2, np.pi/2, np.pi/2)

prm.rho1p = col(0.00622253616772498, 359.132273856264, 0)   # body.CM QP - body.Ai QP (Ai 프레임)
prm.C11 = ang2mat(np.pi, np.pi/2, np.pi/2)

prm.m1 = 18.6342592049469
Ixx = 604104.103452763;   Ixy = 0.126126031622254
Iyy = 565728.29220829;    Iyz = -35.7293476271758
Izz = 40317.2972174389;   Izx = -3.51625272013543e-13
prm.J1p = np.array([[Ixx, Ixy, Izx],
                    [Ixy, Iyy, Iyz],
                    [Izx, Iyz, Izz]])

# ---- system ----
prm.g = -9806.65            # mm/s^2  (RecurDyn KGRAV 와 동일)
prm.F_ex = 10*1000          # cart_pole_16 : FY = step5(time, 0, 10, 1, -10) [N]
                            #   길이가 mm 라 내부 힘 단위는 kg*mm/s^2 -> N 값에 1000 을 곱한다
prm.free = [2, 7]           # 살릴 자유도 : base 전역 Y 병진 + 회전 조인트 (MATLAB 과 같은 1-base)

h = 0.001
t_e = 2

# ======================= initial condition =======================
r0 = col(0, 0, -50)                             # base.Ai 의 전역 위치
p0 = mat2ep(ang2mat(0, -np.pi/2, 0))            # base.Ai 의 전역 자세
q1 = 0.0
dr0 = col(0, 0, 0)
w0 = col(0, 0, 0)
dq1 = 3.0                                       # RecurDyn cart_pole_16 에는 초기속도가 없다

Y = np.block([[r0], [p0], [col(q1)], [dr0], [w0], [col(dq1)]])

# ======================= RK4 =======================
n = round(t_e/h)
T = np.arange(n + 1)*h          # 부동소수 누적 오차 없이 정확한 격자
YY = np.zeros((15, n + 1))
AA = np.zeros((15, n + 1))      # Yp 이력 (가속도 비교용)

YY[:, 0:1] = Y

for k in range(n):              # MATLAB k = 1:n
    t = T[k]
    k1 = dYdt(t,       Y,            prm)
    k2 = dYdt(t + h/2, Y + h/2*k1,   prm)
    k3 = dYdt(t + h/2, Y + h/2*k2,   prm)
    k4 = dYdt(t + h,   Y + h*k3,     prm)

    AA[:, k:k+1] = k1
    Y = Y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

    Y[3:7] = Y[3:7]/np.linalg.norm(Y[3:7])      # 오일러 파라미터 정규화

    YY[:, k+1:k+2] = Y

AA[:, -1:] = dYdt(T[-1], YY[:, -1:], prm)

# ======================= post processing =======================
# cart_pole_16 의 cart_px/py/pz = dx/dy/dz(base.Ai, Ground.origin, Ground.origin)
# Ground.origin 은 QP=(0,0,0), REULER=(0,0,0) 즉 전역계라 변환이 필요없다.
cart_p = YY[0:3, :]         # [cart_px; cart_py; cart_pz]  [mm]
cart_v = YY[8:11, :]        # [cart_vx; cart_vy; cart_vz]  [mm/s]
cart_a = AA[8:11, :]        # [cart_accx; ...]             [mm/s^2]

pend_q = YY[7, :]           # pendulum_q   = az(body.Ai, base.Cij)   [rad]
pend_qd = YY[14, :]         # pendulum_qd                            [rad/s]
pend_qdd = AA[14, :]        # pendulum_qdd                           [rad/s^2]

print(f't = {T[0]:.3f} ~ {T[-1]:.3f} s, {n} steps (h = {h:g})')
print(f'quaternion norm drift : {np.abs(np.linalg.norm(YY[3:7, :], axis=0) - 1).max():.3e}')

# ======================= RecurDyn 비교 =======================
# rec_data.csv : 헤더 없음, 8 열
#   [ index, time, cart_py, cart_vy, cart_accy, pendulum_q, pendulum_qd, pendulum_qdd ]
csv = Path(__file__).resolve().parent / '..' / 'recurdyn' / '01_cart_pole' / 'rec_data.csv'

lab = ['cart_py [mm]',     'cart_vy [mm/s]',      'cart_accy [mm/s^2]',
       'pendulum_q [rad]', 'pendulum_qd [rad/s]', 'pendulum_qdd [rad/s^2]']
mine = np.vstack([cart_p[1, :], cart_v[1, :], cart_a[1, :], pend_q, pend_qd, pend_qdd])

ref = None
if csv.is_file():
    R = np.loadtxt(csv, delimiter=',')
    ref = SimpleNamespace(t=R[:, 1], y=R[:, 2:8].T)
    print(f'\nRecurDyn : rec_data.csv ({ref.t.size} points, t = {ref.t[0]:.3f} ~ {ref.t[-1]:.3f})')
    print(f'{"channel":<24} {"max|err|":>12} {"RMS":>12} {"rel.RMS":>10}')
    err = np.zeros((6, ref.t.size))
    for i in range(6):
        ip = np.interp(ref.t, T, mine[i, :])
        err[i, :] = ip - ref.y[i, :]
        sc = np.abs(ref.y[i, :]).max()
        print(f'{lab[i]:<24} {np.abs(err[i, :]).max():12.4e} {rms(err[i, :]):12.4e} '
              f'{rms(err[i, :])/sc:10.2e}')

    # RecurDyn 의 acc 채널은 알고리즘 감쇠(NDAMPING)로 오염되어 있다.
    # 자기 속도를 미분한 값과 비교해야 실제 모델 일치도가 보인다.
    dt = ref.t[1] - ref.t[0]
    anm = ['cart_accy', 'pendulum_qdd']
    for j, c in enumerate([0, 3]):          # MATLAB 의 cc = [1 4] 를 0-base 로
        a_num = np.gradient(ref.y[c+1, :], dt)          # RecurDyn 의 d(vel)/dt
        am = np.interp(ref.t, T, mine[c+2, :])
        sc = np.abs(ref.y[c+2, :]).max()
        print(f'  {anm[j]:<14} vs 보고 acc {rms(am - ref.y[c+2, :])/sc:8.2e} '
              f'| vs d(vel)/dt {rms(am - a_num)/sc:8.2e}  (rel.RMS)')
else:
    print(f'\n[비교 생략] {csv} 없음')

# ======================= plot =======================
if __name__ == '__main__':
    import os

    import matplotlib

    # WSL/리눅스에서 DISPLAY 가 없으면 창을 띄울 수 없으니 PNG 로 떨군다.
    headless = os.name != 'nt' and not (os.environ.get('DISPLAY') or
                                        os.environ.get('WAYLAND_DISPLAY'))
    if headless:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 2, figsize=(11.5, 7.5), num='MATLAB vs RecurDyn')
    pos = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]   # 좌열 cart, 우열 pendulum
    for i in range(6):
        a = ax[pos[i]]
        a.plot(T, mine[i, :], 'b', linewidth=1.4, label='Python')
        if ref is not None:
            a.plot(ref.t, ref.y[i, :], 'r--', linewidth=1.2, label='RecurDyn')
        a.grid(True); a.set_xlabel('time [s]'); a.set_ylabel(lab[i])
        if i == 0 and ref is not None:
            a.legend(loc='best')
    fig.tight_layout()

    fig2 = None
    if ref is not None:
        fig2, ax2 = plt.subplots(3, 2, figsize=(11.5, 7.5), num='error (Python - RecurDyn)')
        for i in range(6):
            a = ax2[pos[i]]
            a.plot(ref.t, err[i, :], 'k', linewidth=1.1)
            a.grid(True); a.set_xlabel('time [s]'); a.set_ylabel('d ' + lab[i])
        fig2.tight_layout()

    if headless:
        out = Path(__file__).resolve().parent
        fig.savefig(out / 'compare.png', dpi=110)
        print()
        print(f'[headless] {out / "compare.png"} 저장')
        if fig2 is not None:
            fig2.savefig(out / 'error.png', dpi=110)
            print(f'[headless] {out / "error.png"} 저장')
    else:
        plt.show()
