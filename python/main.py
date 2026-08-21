r"""main.m 1:1 변환.

실행:  .venv\Scripts\python.exe main.py
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from dYdt import dYdt
from model import create_parameters, full_state
from util import rms

# ======================= parameter =======================
# 단위계 : SI (m, kg, s, N). 관성모멘트 단위는 kg*m^2 이다.
prm = create_parameters()

h = 0.001
t_e = 2

# ======================= initial condition =======================
Y = full_state([0.0, 0.0, 0.0, 3.0])  # dq1=3은 RecurDyn의 초기속도 0과 다르다.

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
cart_p = YY[0:3, :]         # [cart_px; cart_py; cart_pz]  [m]
cart_v = YY[8:11, :]        # [cart_vx; cart_vy; cart_vz]  [m/s]
cart_a = AA[8:11, :]        # [cart_accx; ...]             [m/s^2]

pend_q = YY[7, :]           # pendulum_q   = az(body.Ai, base.Cij)   [rad]
pend_qd = YY[14, :]         # pendulum_qd                            [rad/s]
pend_qdd = AA[14, :]        # pendulum_qdd                           [rad/s^2]

print(f't = {T[0]:.3f} ~ {T[-1]:.3f} s, {n} steps (h = {h:g})')
print(f'quaternion norm drift : {np.abs(np.linalg.norm(YY[3:7, :], axis=0) - 1).max():.3e}')

# ======================= RecurDyn 비교 =======================
# rec_data.csv : 헤더 없음, 8 열 (RecurDyn 원본 병진 채널은 mm 단위)
#   [ index, time, cart_py, cart_vy, cart_accy, pendulum_q, pendulum_qd, pendulum_qdd ]
csv = Path(__file__).resolve().parent / '..' / 'recurdyn' / '01_cart_pole' / 'rec_data.csv'

lab = ['cart_py [m]',      'cart_vy [m/s]',       'cart_accy [m/s^2]',
       'pendulum_q [rad]', 'pendulum_qd [rad/s]', 'pendulum_qdd [rad/s^2]']
mine = np.vstack([cart_p[1, :], cart_v[1, :], cart_a[1, :], pend_q, pend_qd, pend_qdd])

ref = None
if csv.is_file():
    R = np.loadtxt(csv, delimiter=',')
    ref = SimpleNamespace(t=R[:, 1], y=R[:, 2:8].T)
    ref.y[0:3, :] *= 1e-3                       # mm 계열 -> SI(m 계열)
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

    # 리눅스에서 디스플레이가 없으면 창을 띄울 수 없으므로 Agg 로 내린다.
    if os.name != 'nt' and not (os.environ.get('DISPLAY') or
                                os.environ.get('WAYLAND_DISPLAY')):
        matplotlib.use('Agg')

    # 최종적으로 잡힌 backend 가 비대화형이면 PNG 로 떨군다.
    # (MPLBACKEND=Agg 로 강제했거나 GUI 툴킷이 없어 Agg 로 떨어진 경우까지 포함)
    headless = matplotlib.get_backend().lower() in {
        'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template'}
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
