"""MATLAB 에 내장되어 있어 .m 에는 대응 파일이 없는 것들만 모아둔 보조 모듈."""

import numpy as np


def col(*vals):
    """MATLAB 의 [a; b; c] 열벡터.  col(0, 0, -50) -> shape (3, 1)"""
    return np.array(vals, dtype=float).reshape(-1, 1)


def rms(x):
    """MATLAB 의 rms()"""
    x = np.asarray(x, dtype=float).ravel()
    return float(np.sqrt(np.mean(x**2)))
