"""Cart-pendulum 모델의 SI 파라미터와 상태 변환."""

from types import SimpleNamespace

import numpy as np

from ang2mat import ang2mat
from mat2ep import mat2ep
from util import col


BASE_Z = -0.05
BASE_ATTITUDE = mat2ep(ang2mat(0.0, -np.pi/2, 0.0))


def create_parameters():
    """dYdt에서 사용하는 SI 단위계 파라미터를 생성한다."""
    prm = SimpleNamespace()

    # base body (기준계 = base.Ai)
    prm.rho0p = col(0.0, 0.0, 0.0)
    prm.C00 = ang2mat(0.0, np.pi/2, 0.0)
    prm.m0 = 111.015764646288
    prm.J0p = np.array([
        [9.34382685772924, 0.0, 0.0],
        [0.0, 0.27753941161572, 0.0],
        [0.0, 0.0, 9.43633999493448],
    ])

    # link body (기준계 = body.Ai)
    prm.s01p = col(0.0, 0.0, 0.0)
    prm.C01 = ang2mat(np.pi/2, np.pi/2, np.pi/2)
    prm.rho1p = col(6.22253616772498e-6, 0.359132273856264, 0.0)
    prm.C11 = ang2mat(np.pi, np.pi/2, np.pi/2)
    prm.m1 = 18.6342592049469
    prm.J1p = np.array([
        [0.604104103452763, 1.26126031622254e-7, -3.51625272013543e-19],
        [1.26126031622254e-7, 0.56572829220829, -3.57293476271758e-5],
        [-3.51625272013543e-19, -3.57293476271758e-5, 0.0403172972174389],
    ])

    # system
    prm.g = -9.80665
    prm.F_ex = 10.0
    prm.free = [2, 7]  # base 전역 Y 병진 + 회전 조인트 (MATLAB과 같은 1-base)
    return prm


def full_state(x):
    """축약 절대상태 [y, q1, dy, dq1]를 dYdt의 15차원 상태로 변환한다."""
    y, q1, dy, dq1 = np.asarray(x, dtype=float).reshape(4)
    r0 = col(0.0, y, BASE_Z)
    dr0 = col(0.0, dy, 0.0)
    w0 = col(0.0, 0.0, 0.0)
    return np.block([
        [r0],
        [BASE_ATTITUDE.copy()],
        [col(q1)],
        [dr0],
        [w0],
        [col(dq1)],
    ])


def reduced_state(Y):
    """dYdt의 15차원 상태에서 [y, q1, dy, dq1]를 추출한다."""
    Y = np.asarray(Y, dtype=float).reshape(15, -1)
    return np.array([Y[1, 0], Y[7, 0], Y[9, 0], Y[14, 0]])
