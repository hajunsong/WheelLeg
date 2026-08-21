"""dYdt.m 1:1 변환."""

import numpy as np

from step5 import step5
from tilde import tilde
from util import col


def dYdt(t, Y, prm, full=False, Fy=None):
    """
    상태벡터 Y 를 받아 미분 Yp 를 돌려준다. (RK4 의 각 stage 에서 호출)

      Y  = [ r0(3) ; p0(4) ; q1 ;  dr0(3) ;  w0(3) ; dq1  ]   (15x1)
      Yp = [ dr0(3); dp0(4); dq1;  ddr0(3);  dw0(3); ddq1 ]   (15x1)

    full=True 이면 (Yp, out) 을 돌려준다.  (MATLAB 의 두번째 출력)
    Fy가 주어지면 시간 함수 대신 해당 제어력[N]을 사용한다.
    """
    # ---------------- Y2qdq ----------------
    # MATLAB 은 1-base, Python 은 0-base 라 인덱스만 한 칸씩 당겨진다.
    r0 = Y[0:3]                 # MATLAB Y(1:3)
    p0 = Y[3:7]                 # MATLAB Y(4:7)
    q1 = Y[7, 0]                # MATLAB Y(8)
    dr0 = Y[8:11]               # MATLAB Y(9:11)
    w0 = Y[11:14]               # MATLAB Y(12:14)
    dq1 = Y[14, 0]              # MATLAB Y(15)

    # RK4 중간 stage 에서 노름이 흐르므로 매번 투영해준다
    p0 = p0 / np.linalg.norm(p0)

    # ---------------- base body ----------------
    e0 = p0[1:4]
    E0 = np.block([-e0,  tilde(e0) + p0[0, 0]*np.eye(3)])
    G0 = np.block([-e0, -tilde(e0) + p0[0, 0]*np.eye(3)])
    A0 = E0 @ G0.T

    J0c = A0 @ prm.C00 @ prm.J0p @ (A0 @ prm.C00).T

    rho0 = A0 @ prm.rho0p
    r0c = r0 + rho0

    w0t = tilde(w0)
    r0t = tilde(r0)

    dr0c = dr0 + w0t @ rho0

    dr0t = tilde(dr0)
    dr0ct = tilde(dr0c)
    r0ct = tilde(r0c)

    Y0h = np.block([[dr0 + r0t @ w0], [w0]])

    # 외력 : RecurDyn TRANSLATIONAL_FORCE, FY = step5(time, 0, F, 1, -F) [N], RM = global
    #        작용점 base.Marker4 가 base.CM 과 같은 위치라 순수 CM 힘으로 들어간다
    if Fy is None:
        Fy = step5(t, 0.0, prm.F_ex, 1.0, -prm.F_ex)
    else:
        Fy = float(Fy)

    f0c = col(0.0, Fy, prm.m0*prm.g)
    t0c = col(0.0, 0.0, 0.0)

    M0h = np.block([[prm.m0*np.eye(3), -prm.m0*r0ct],
                    [prm.m0*r0ct,       J0c - prm.m0*r0ct @ r0ct]])
    Q0h = np.block([[f0c + prm.m0*dr0ct @ w0],
                    [t0c + r0ct @ f0c + prm.m0*r0ct @ dr0ct @ w0 - w0t @ J0c @ w0]])

    # ---------------- link body ----------------
    A01pp = np.array([[np.cos(q1), -np.sin(q1), 0.0],
                      [np.sin(q1),  np.cos(q1), 0.0],
                      [       0.0,         0.0, 1.0]])
    A1 = A0 @ prm.C01 @ A01pp
    s01 = A1 @ prm.s01p
    r1 = r0 + s01

    J1c = A1 @ prm.C11 @ prm.J1p @ (A1 @ prm.C11).T

    rho1 = A1 @ prm.rho1p
    r1c = r1 + rho1
    r1ct = tilde(r1c)

    H1 = A0 @ prm.C01 @ col(0.0, 0.0, 1.0)
    w1 = w0 + H1*dq1
    w1t = tilde(w1)
    dr1 = dr0 + w0t @ s01
    r1t = tilde(r1)

    B1 = np.block([[r1t @ H1], [H1]])

    dr1t = tilde(dr1)
    dr1c = dr1 + w1t @ rho1
    dr1ct = tilde(dr1c)

    dH1 = w0t @ H1
    D1 = np.block([[dr1t @ H1 + r1t @ dH1], [dH1]])*dq1
    Y1h = Y0h + B1*dq1

    f1c = col(0.0, 0.0, prm.m1*prm.g)
    t1c = col(0.0, 0.0, 0.0)

    M1h = np.block([[prm.m1*np.eye(3), -prm.m1*r1ct],
                    [prm.m1*r1ct,       J1c - prm.m1*r1ct @ r1ct]])
    Q1h = np.block([[f1c + prm.m1*dr1ct @ w1],
                    [t1c + r1ct @ f1c + prm.m1*r1ct @ dr1ct @ w1 - w1t @ J1c @ w1]])

    # ---------------- mass / force ----------------
    K1 = M1h
    K0 = K1 + M0h

    L1 = Q1h
    L0 = L1 + Q0h - K1 @ D1

    # ---------------- EQM ----------------
    M = np.block([[K0,       K1 @ B1],
                  [B1.T @ K1, B1.T @ K1 @ B1]])
    Q = np.block([[L0],
                  [B1.T @ (L1 - K1 @ D1)]])

    # 구속 : prm.free 이외의 자유도는 가속도 0 으로 잠근다
    # prm.free 는 MATLAB 과 같은 1-base 로 두고 여기서만 0-base 로 바꾼다.
    lock = np.array([i for i in range(7) if (i + 1) not in prm.free])
    M[lock, :] = 0.0
    M[lock, lock] = 1.0             # MATLAB sub2ind 대각 대입과 동일
    Q[lock] = 0.0

    ddq = np.linalg.solve(M, Q)
    dY0h = ddq[0:6]
    ddq1 = ddq[6, 0]

    # ---------------- base body acceleration ----------------
    dp0 = 0.5*E0.T @ w0

    T0 = np.block([[np.eye(3),   -r0t],
                   [np.zeros((3, 3)), np.eye(3)]])
    R0 = np.block([[dr0t @ w0], [np.zeros((3, 1))]])
    dY0b = T0 @ dY0h - R0

    ddr0 = dY0b[0:3]
    dw0 = dY0b[3:6]

    dw0t = tilde(dw0)
    ddr0c = ddr0 + dw0t @ rho0 + w0t @ w0t @ rho0

    # ---------------- link body acceleration ----------------
    dY1h = dY0h + B1*ddq1 + D1

    T1 = np.block([[np.eye(3),   -r1t],
                   [np.zeros((3, 3)), np.eye(3)]])
    dT1 = np.block([[np.zeros((3, 3)), -dr1t],
                    [np.zeros((3, 3)),  np.zeros((3, 3))]])
    dY1b = dT1 @ Y1h + T1 @ dY1h

    ddr1 = dY1b[0:3]
    dw1 = dY1b[3:6]
    dw1t = tilde(dw1)

    ddr1c = ddr1 + dw1t @ rho1 + w1t @ w1t @ rho1

    # ---------------- dqddq2Yp ----------------
    Yp = np.block([[dr0], [dp0], [col(dq1)], [ddr0], [dw0], [col(ddq1)]])

    if not full:
        return Yp

    out = dict(A0=A0, A1=A1, r0c=r0c, r1c=r1c,
               dr0c=dr0c, dr1c=dr1c, ddr0c=ddr0c, ddr1c=ddr1c,
               w1=w1, dw0=dw0, dw1=dw1, ddq1=ddq1,
               Fy=Fy, M=M, Q=Q)
    return Yp, out
