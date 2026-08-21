"""비선형 cart-pendulum 모델의 직립 평형점 LQR 제어.

순서:
  비선형 EOM -> 평형점 -> 수치 선형화 -> 제어 가능성 확인
  -> Q/R 선정 -> CARE 풀이 -> 폐루프 극 확인 -> 비선형 폐루프 시뮬레이션
"""

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_are

from dYdt import dYdt
from model import create_parameters, full_state


def nonlinear_dynamics(t, x, prm, force):
    """축약 절대상태 x=[y, q1, dy, dq1]의 비선형 상태방정식."""
    Yp = dYdt(t, full_state(x), prm, Fy=force)
    return np.array([Yp[1, 0], Yp[7, 0], Yp[9, 0], Yp[14, 0]])


def upright_equilibrium(prm, y_ref=0.0):
    """링크 CM이 조인트 바로 위에 놓이는 직립 평형점 (xe, ue)."""
    a = float(prm.rho1p[0, 0])
    b = float(prm.rho1p[1, 0])
    q1e = np.arctan2(a, b)
    return np.array([y_ref, q1e, 0.0, 0.0]), 0.0


def numerical_linearization(fun, xe, ue):
    """중앙차분으로 A=df/dx, B=df/du를 계산한다."""
    n = xe.size
    A = np.zeros((n, n))
    dx = np.array([1e-6, 1e-7, 1e-6, 1e-7])
    for i in range(n):
        xp = xe.copy()
        xm = xe.copy()
        xp[i] += dx[i]
        xm[i] -= dx[i]
        A[:, i] = (fun(xp, ue) - fun(xm, ue))/(2.0*dx[i])

    du = 1e-4
    B = ((fun(xe, ue + du) - fun(xe, ue - du))/(2.0*du)).reshape(n, 1)
    return A, B


def controllability_matrix(A, B):
    """C=[B AB ... A^(n-1)B]."""
    blocks = [B]
    for _ in range(1, A.shape[0]):
        blocks.append(A @ blocks[-1])
    return np.hstack(blocks)


def is_stabilizable(A, B, tol=1e-9):
    """불안정/중립 고유모드에 PBH rank 조건을 적용한다."""
    n = A.shape[0]
    for pole in np.linalg.eigvals(A):
        if pole.real >= -tol:
            pbh = np.hstack([pole*np.eye(n) - A, B])
            if np.linalg.matrix_rank(pbh, tol=tol) < n:
                return False
    return True


def lqr_gain(A, B, Q, R):
    """연속시간 CARE의 안정화 해 P와 LQR 이득 K를 계산한다."""
    P = solve_continuous_are(A, B, Q, R)
    K = np.linalg.solve(R, B.T @ P)
    return K, P


def angle_error(q, q_ref):
    """각도 오차를 [-pi, pi] 범위로 감싼다."""
    return np.arctan2(np.sin(q - q_ref), np.cos(q - q_ref))


def main():
    prm = create_parameters()
    y_ref0 = 0.0
    v_ref = -1.0                 # cart 목표 속도 [m/s]
    xe, ue = upright_equilibrium(prm, y_ref=y_ref0)

    def reference(t):
        """등속 이동 기준궤적 [y_ref(t), q1e, v_ref, 0]."""
        return np.array([y_ref0 + v_ref*t, xe[1], v_ref, 0.0])

    def f(x, u):
        return nonlinear_dynamics(0.0, x, prm, u)

    equilibrium_residual = f(xe, ue)
    A, B = numerical_linearization(f, xe, ue)

    C = controllability_matrix(A, B)
    controllability_rank = np.linalg.matrix_rank(C)
    stabilizable = is_stabilizable(A, B)
    if controllability_rank != A.shape[0]:
        raise RuntimeError(f'선형화 모델이 제어 가능하지 않습니다: rank={controllability_rank}')
    if not stabilizable:
        raise RuntimeError('선형화 모델이 stabilizable하지 않습니다.')

    # Bryson 규칙: 각 상태/입력의 허용 최대값 역제곱을 가중치로 사용한다.
    y_max = 0.20                  # m
    q_max = np.deg2rad(5.0)      # rad
    dy_max = 1.0                 # m/s
    dq_max = 2.0                 # rad/s
    force_scale = 100.0          # N, 비용함수 설계 기준
    force_limit = 500.0          # N, 실제 액추에이터 포화
    Q = np.diag(1.0/np.square([y_max, q_max, dy_max, dq_max]))
    R = np.array([[1.0/force_scale**2]])

    K, P = lqr_gain(A, B, Q, R)
    closed_loop_poles = np.linalg.eigvals(A - B @ K)

    # 그림의 Hamiltonian 절차는 CARE 해법의 이론적 배경이다.
    # 실제 계산은 고유벡터 분할보다 수치적으로 안정적인 solve_continuous_are를 사용한다.
    Rinv = np.linalg.inv(R)
    H = np.block([[A, -B @ Rinv @ B.T],
                  [-Q, -A.T]])
    hamiltonian_poles = np.linalg.eigvals(H)
    stable_hamiltonian_count = np.count_nonzero(hamiltonian_poles.real < 0.0)

    print('선형화 직립 평형점 xe =', xe)
    print(f'평형 입력 ue = {ue:.6g} N')
    print('평형점 잔차 f(xe,ue) =', equilibrium_residual)
    print(f'cart 기준궤적: y_ref(t) = {y_ref0:g} + ({v_ref:g})*t [m]')
    print('\nA =\n', A)
    print('\nB =\n', B)
    print(f'\ncontrollability rank = {controllability_rank}/{A.shape[0]}')
    print('stabilizable =', stabilizable)
    print('\nQ =\n', Q)
    print('\nR =\n', R)
    print('\nK =\n', K)
    print('\neig(A-BK) =\n', closed_loop_poles)
    print(f'Hamiltonian 안정 고유값 개수 = {stable_hamiltonian_count}/{A.shape[0]}')

    def feedback_force(t, x):
        x_ref = reference(t)
        error = x - x_ref
        error[1] = angle_error(x[1], x_ref[1])
        command = ue - (K @ error).item()
        return np.clip(command, -force_limit, force_limit)

    def closed_loop(t, x):
        return nonlinear_dynamics(t, x, prm, feedback_force(t, x))

    # LQR은 국소 제어기이므로 직립점 부근(5 deg)에서 시작한다.
    x0 = xe + np.array([0.00, np.deg2rad(0.0), 0.0, 0.5])
    t_eval = np.linspace(0.0, 10.0, 2001)
    sol = solve_ivp(
        closed_loop, (t_eval[0], t_eval[-1]), x0, t_eval=t_eval,
        method='RK45', rtol=1e-8, atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    x_ref = np.column_stack([reference(t) for t in sol.t])
    force = np.array([feedback_force(t, x) for t, x in zip(sol.t, sol.y.T)])
    final_error = sol.y[:, -1] - x_ref[:, -1]
    final_error[1] = angle_error(sol.y[1, -1], x_ref[1, -1])
    print('\n최종 추종오차 =', final_error)
    print(f'최대 제어력 = {np.max(np.abs(force)):.3f} N')
    animate_response(sol.t, sol.y, x_ref, force, prm)
    plot_response(sol.t, sol.y, x_ref, force, closed_loop_poles)


def animate_response(t, x, x_ref, force, prm):
    """카트와 진자의 비선형 폐루프 응답을 GIF로 저장한다."""
    import os
    import matplotlib

    if os.name != 'nt' and not (
            os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.patches import Circle, Rectangle

    fps = 30
    frame_count = min(int(np.ceil((t[-1] - t[0])*fps)) + 1, t.size)
    frame_indices = np.linspace(0, t.size - 1, frame_count, dtype=int)

    length = float(np.linalg.norm(prm.rho1p))
    theta = x[1] - x_ref[1]
    bob_y = x[0] - length*np.sin(theta)
    bob_z = length*np.cos(theta)

    cart_width = 0.24
    cart_height = 0.10
    view_half_width = max(0.6, 1.5*(length + cart_width))

    fig, ax = plt.subplots(figsize=(9.0, 4.8), num='Nonlinear LQR animation')
    ax.set_xlim(x_ref[0, 0] - view_half_width,
                x_ref[0, 0] + view_half_width)
    ax.set_ylim(-0.18, length + 0.18)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('global Y [m]')
    ax.set_ylabel('global Z [m]')
    ax.set_title('Cart-pendulum nonlinear LQR control')
    ax.grid(True, alpha=0.25)
    ax.axhline(-cart_height, color='0.25', linewidth=2.0)
    target_line = ax.axvline(
        x_ref[0, 0], color='0.5', linestyle='--', linewidth=1.0,
    )

    cart = Rectangle((0.0, -cart_height), cart_width, cart_height,
                     facecolor='#377eb8', edgecolor='black', zorder=3)
    bob = Circle((0.0, 0.0), 0.025, facecolor='#e41a1c',
                 edgecolor='black', zorder=5)
    ax.add_patch(cart)
    ax.add_patch(bob)
    rod, = ax.plot([], [], color='black', linewidth=3.0, zorder=4)
    trail, = ax.plot([], [], color='#e41a1c', linewidth=1.0, alpha=0.35)
    force_line, = ax.plot([], [], color='#4daf4a', linewidth=4.0)
    info = ax.text(0.02, 0.96, '', transform=ax.transAxes, va='top',
                   family='monospace')

    force_scale = max(np.max(np.abs(force)), 1.0)

    def update(frame):
        i = frame_indices[frame]
        cart_y = x[0, i]
        target_y = x_ref[0, i]
        ax.set_xlim(target_y - view_half_width, target_y + view_half_width)
        target_line.set_xdata([target_y, target_y])
        cart.set_xy((cart_y - cart_width/2.0, -cart_height))
        rod.set_data([cart_y, bob_y[i]], [0.0, bob_z[i]])
        bob.center = (bob_y[i], bob_z[i])

        start = max(0, i - round(t.size/(t[-1] - t[0])))
        trail.set_data(bob_y[start:i+1], bob_z[start:i+1])

        arrow_length = 0.18*force[i]/force_scale
        force_line.set_data([cart_y, cart_y + arrow_length],
                            [-0.5*cart_height, -0.5*cart_height])
        info.set_text(
            f't     = {t[i]:5.2f} s\n'
            f'y     = {cart_y:+7.3f} m\n'
            f'y_ref = {target_y:+7.3f} m\n'
            f'angle = {np.rad2deg(theta[i]):+7.3f} deg\n'
            f'force = {force[i]:+7.1f} N'
        )
        return cart, rod, bob, trail, force_line, target_line, info

    animation = FuncAnimation(
        fig, update, frames=frame_count, interval=1000/fps, blit=False,
    )
    output = Path(__file__).resolve().parent / 'lqr_animation.gif'
    animation.save(output, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    print(f'[animation] {output} 저장')


def plot_response(t, x, x_ref, force, closed_loop_poles):
    """폐루프 응답을 화면에 표시하거나 lqr_response.png로 저장한다."""
    import os
    import matplotlib

    if os.name != 'nt' and not (
            os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 2, figsize=(11.5, 8.0), num='Nonlinear LQR response')
    ax[0, 0].plot(t, x[0], linewidth=1.4, label='actual')
    ax[0, 0].plot(t, x_ref[0], 'k--', linewidth=0.9, label='reference')
    ax[0, 0].set_ylabel('cart y [m]')
    ax[0, 0].legend(loc='best')

    angle_err = np.arctan2(
        np.sin(x[1] - x_ref[1]), np.cos(x[1] - x_ref[1]),
    )
    ax[0, 1].plot(t, np.rad2deg(angle_err), linewidth=1.4)
    ax[0, 1].axhline(0.0, color='k', linestyle='--', linewidth=0.8)
    ax[0, 1].set_ylabel('pendulum error [deg]')

    ax[1, 0].plot(t, x[2], linewidth=1.4, label='actual')
    ax[1, 0].plot(t, x_ref[2], 'k--', linewidth=0.9, label='reference')
    ax[1, 0].set_ylabel('cart velocity [m/s]')
    ax[1, 0].legend(loc='best')

    ax[1, 1].plot(t, x[3], linewidth=1.4)
    ax[1, 1].plot(t, x_ref[3], 'k--', linewidth=0.9)
    ax[1, 1].set_ylabel('pendulum rate [rad/s]')

    ax[2, 0].plot(t, force, linewidth=1.4)
    ax[2, 0].set_ylabel('control force [N]')

    ax[2, 1].axis('off')
    poles = '\n'.join(f'{p.real:+.3f}{p.imag:+.3f}j' for p in closed_loop_poles)
    ax[2, 1].text(0.05, 0.95, 'eig(A-BK)\n' + poles, va='top', family='monospace')

    for a in ax.flat:
        if a.axison:
            a.grid(True)
            a.set_xlabel('time [s]')
    fig.tight_layout()

    if matplotlib.get_backend().lower() in {
            'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template'}:
        output = Path(__file__).resolve().parent / 'lqr_response.png'
        fig.savefig(output, dpi=120)
        print(f'\n[headless] {output} 저장')
    else:
        plt.show()


if __name__ == '__main__':
    main()
