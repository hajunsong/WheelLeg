"""tilde.m 1:1 변환"""

import numpy as np


def tilde(x):
    x = np.asarray(x, dtype=float).ravel()
    return np.array([[    0.0, -x[2],  x[1]],
                     [   x[2],   0.0, -x[0]],
                     [  -x[1],  x[0],   0.0]])
