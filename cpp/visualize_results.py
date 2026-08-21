"""C++ CSV 결과를 PNG와 GIF로 시각화한다."""

import argparse
import os
from pathlib import Path

import matplotlib
import numpy as np

if os.name != 'nt' and not (
        os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, Rectangle


HERE = Path(__file__).resolve().parent


def plot_open_loop(data):
    labels = [
        ('cart_y', 'cart y [m]'),
        ('q1', 'pendulum q1 [rad]'),
        ('cart_vy', 'cart velocity [m/s]'),
        ('dq1', 'pendulum rate [rad/s]'),
        ('cart_ay', 'cart acceleration [m/s²]'),
        ('ddq1', 'pendulum acceleration [rad/s²]'),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(11.5, 8.0), num='C++ open-loop')
    for axis, (field, label) in zip(axes.flat, labels):
        axis.plot(data['time'], data[field], linewidth=1.3)
        axis.set_xlabel('time [s]')
        axis.set_ylabel(label)
        axis.grid(True)
    fig.tight_layout()
    output = HERE / 'cpp_open_loop.png'
    fig.savefig(output, dpi=120)
    plt.close(fig)
    print(f'[plot] {output} 저장')


def plot_lqr(data):
    t = data['time']
    angle_error = np.arctan2(
        np.sin(data['q1'] - data['q1_ref']),
        np.cos(data['q1'] - data['q1_ref']),
    )
    fig, ax = plt.subplots(3, 2, figsize=(11.5, 8.0), num='C++ nonlinear LQR')

    ax[0, 0].plot(t, data['y'], label='actual')
    ax[0, 0].plot(t, data['y_ref'], 'k--', linewidth=0.9, label='reference')
    ax[0, 0].set_ylabel('cart y [m]')
    ax[0, 0].legend(loc='best')

    ax[0, 1].plot(t, np.rad2deg(angle_error))
    ax[0, 1].axhline(0.0, color='k', linestyle='--', linewidth=0.8)
    ax[0, 1].set_ylabel('pendulum error [deg]')

    ax[1, 0].plot(t, data['dy'], label='actual')
    ax[1, 0].plot(t, data['dy_ref'], 'k--', linewidth=0.9, label='reference')
    ax[1, 0].set_ylabel('cart velocity [m/s]')
    ax[1, 0].legend(loc='best')

    ax[1, 1].plot(t, data['dq1'])
    ax[1, 1].plot(t, data['dq1_ref'], 'k--', linewidth=0.9)
    ax[1, 1].set_ylabel('pendulum rate [rad/s]')

    ax[2, 0].plot(t, data['force'])
    ax[2, 0].set_ylabel('control force [N]')
    ax[2, 1].axis('off')

    for axis in ax.flat:
        if axis.axison:
            axis.set_xlabel('time [s]')
            axis.grid(True)
    fig.tight_layout()
    output = HERE / 'cpp_lqr_response.png'
    fig.savefig(output, dpi=120)
    plt.close(fig)
    print(f'[plot] {output} 저장')


def animate_lqr(data):
    t = data['time']
    fps = 30
    frame_count = min(int(np.ceil((t[-1] - t[0])*fps)) + 1, t.size)
    indices = np.linspace(0, t.size - 1, frame_count, dtype=int)

    length = np.hypot(6.22253616772498e-6, 0.359132273856264)
    theta = data['q1'] - data['q1_ref']
    bob_y = data['y'] - length*np.sin(theta)
    bob_z = length*np.cos(theta)
    cart_width = 0.24
    cart_height = 0.10
    view_half_width = max(0.6, 1.5*(length + cart_width))

    fig, ax = plt.subplots(figsize=(9.0, 4.8), num='C++ nonlinear LQR animation')
    ax.set_xlim(data['y_ref'][0] - view_half_width,
                data['y_ref'][0] + view_half_width)
    ax.set_ylim(-0.18, length + 0.18)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('global Y [m]')
    ax.set_ylabel('global Z [m]')
    ax.set_title('C++/Eigen cart-pendulum LQR control')
    ax.grid(True, alpha=0.25)
    ax.axhline(-cart_height, color='0.25', linewidth=2.0)
    target_line = ax.axvline(
        data['y_ref'][0], color='0.5', linestyle='--', linewidth=1.0,
    )

    cart = Rectangle(
        (0.0, -cart_height), cart_width, cart_height,
        facecolor='#377eb8', edgecolor='black', zorder=3,
    )
    bob = Circle(
        (0.0, 0.0), 0.025, facecolor='#e41a1c',
        edgecolor='black', zorder=5,
    )
    ax.add_patch(cart)
    ax.add_patch(bob)
    rod, = ax.plot([], [], color='black', linewidth=3.0, zorder=4)
    trail, = ax.plot([], [], color='#e41a1c', linewidth=1.0, alpha=0.35)
    force_line, = ax.plot([], [], color='#4daf4a', linewidth=4.0)
    info = ax.text(0.02, 0.96, '', transform=ax.transAxes,
                   va='top', family='monospace')
    force_scale = max(np.max(np.abs(data['force'])), 1.0)
    samples_per_second = round(data.size/(t[-1] - t[0]))

    def update(frame):
        i = indices[frame]
        cart_y = data['y'][i]
        target_y = data['y_ref'][i]
        ax.set_xlim(target_y - view_half_width, target_y + view_half_width)
        target_line.set_xdata([target_y, target_y])
        cart.set_xy((cart_y - cart_width/2.0, -cart_height))
        rod.set_data([cart_y, bob_y[i]], [0.0, bob_z[i]])
        bob.center = (bob_y[i], bob_z[i])

        start = max(0, i - samples_per_second)
        trail.set_data(bob_y[start:i+1], bob_z[start:i+1])
        arrow_length = 0.18*data['force'][i]/force_scale
        force_line.set_data(
            [cart_y, cart_y + arrow_length],
            [-0.5*cart_height, -0.5*cart_height],
        )
        info.set_text(
            f't     = {t[i]:5.2f} s\n'
            f'y     = {cart_y:+7.3f} m\n'
            f'y_ref = {target_y:+7.3f} m\n'
            f'angle = {np.rad2deg(theta[i]):+7.3f} deg\n'
            f'force = {data["force"][i]:+7.1f} N'
        )

    animation = FuncAnimation(
        fig, update, frames=frame_count, interval=1000/fps, blit=False,
    )
    output = HERE / 'cpp_lqr_animation.gif'
    animation.save(output, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    print(f'[animation] {output} 저장')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-gif', action='store_true',
                        help='느린 GIF 생성을 생략한다.')
    args = parser.parse_args()

    open_loop_path = HERE / 'open_loop.csv'
    lqr_path = HERE / 'lqr_response.csv'
    if not open_loop_path.is_file() or not lqr_path.is_file():
        raise FileNotFoundError(
            '먼저 wheelleg_open_loop과 wheelleg_lqr을 실행하세요.'
        )

    open_loop = np.genfromtxt(open_loop_path, delimiter=',', names=True)
    lqr = np.genfromtxt(lqr_path, delimiter=',', names=True)
    plot_open_loop(open_loop)
    plot_lqr(lqr)
    if not args.skip_gif:
        animate_lqr(lqr)


if __name__ == '__main__':
    main()
