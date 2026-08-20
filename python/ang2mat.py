"""ang2mat.m 1:1 변환.  RecurDyn REULER (Body 313, ZXZ) -> 회전행렬"""

import numpy as np


def ang2mat(psi, theta, phi):
    Rz1 = np.array([[np.cos(psi), -np.sin(psi), 0.0],
                    [np.sin(psi),  np.cos(psi), 0.0],
                    [        0.0,          0.0, 1.0]])
    Rx = np.array([[1.0,            0.0,             0.0],
                   [0.0, np.cos(theta), -np.sin(theta)],
                   [0.0, np.sin(theta),  np.cos(theta)]])
    Rz2 = np.array([[np.cos(phi), -np.sin(phi), 0.0],
                    [np.sin(phi),  np.cos(phi), 0.0],
                    [        0.0,          0.0, 1.0]])

    mat = Rz1 @ Rx @ Rz2
    return mat
