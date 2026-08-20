"""ep2mat.m 1:1 변환.  오일러 파라미터 [e0;e1;e2;e3] -> 회전행렬.  mat2ep 의 역변환."""

import numpy as np

from tilde import tilde


def ep2mat(p):
    p = np.asarray(p, dtype=float).reshape(-1, 1)
    p = p / np.linalg.norm(p)
    e0 = p[0, 0]
    e = p[1:4]

    E = np.block([-e,  tilde(e) + e0*np.eye(3)])
    G = np.block([-e, -tilde(e) + e0*np.eye(3)])
    A = E @ G.T
    return A
