"""mat2ep.m 1:1 변환.  회전행렬 -> 오일러 파라미터 (단위 쿼터니언, scalar first)"""

import numpy as np


def mat2ep(A, p_prev=None):
    """
    p = mat2ep(A)          A(3x3) -> p = [e0;e1;e2;e3],  e0 >= 0 으로 정규화
    p = mat2ep(A, p_prev)  p_prev 와 부호가 연속이 되도록 선택 (적분 루프용)

    규약 (dYdt 의 E, G 와 동일):
        E = [-e,  tilde(e) + e0*eye(3)]
        G = [-e, -tilde(e) + e0*eye(3)]
        A = E @ G.T

    Shepperd 방법: 4*ei^2 후보 중 가장 큰 것부터 구해 0 나눗셈을 피한다.
    """
    A = np.asarray(A, dtype=float)
    tr = A[0, 0] + A[1, 1] + A[2, 2]

    # 4*e0^2, 4*e1^2, 4*e2^2, 4*e3^2
    c = np.array([1.0 + tr,
                  1.0 + 2.0*A[0, 0] - tr,
                  1.0 + 2.0*A[1, 1] - tr,
                  1.0 + 2.0*A[2, 2] - tr])

    k = int(np.argmax(c))       # MATLAB [cmax, k] = max(c) 와 동일 (첫 최대값)
    s = np.sqrt(c[k])           # 2*|e_k|
    d = 0.5 / s                 # 1/(4*e_k)

    if k == 0:                  # MATLAB case 1
        e0 = 0.5*s
        e1 = (A[2, 1] - A[1, 2])*d
        e2 = (A[0, 2] - A[2, 0])*d
        e3 = (A[1, 0] - A[0, 1])*d
    elif k == 1:                # MATLAB case 2
        e1 = 0.5*s
        e0 = (A[2, 1] - A[1, 2])*d
        e2 = (A[0, 1] + A[1, 0])*d
        e3 = (A[0, 2] + A[2, 0])*d
    elif k == 2:                # MATLAB case 3
        e2 = 0.5*s
        e0 = (A[0, 2] - A[2, 0])*d
        e1 = (A[0, 1] + A[1, 0])*d
        e3 = (A[1, 2] + A[2, 1])*d
    else:                       # MATLAB otherwise
        e3 = 0.5*s
        e0 = (A[1, 0] - A[0, 1])*d
        e1 = (A[0, 2] + A[2, 0])*d
        e2 = (A[1, 2] + A[2, 1])*d

    p = np.array([[e0], [e1], [e2], [e3]])
    p = p / np.linalg.norm(p)

    # p 와 -p 는 같은 자세이므로 부호를 하나 골라야 한다
    if p_prev is None:
        if p[0, 0] < 0.0:
            p = -p                              # e0 >= 0  (theta in [0, 180deg])
    else:
        p_prev = np.asarray(p_prev, dtype=float).reshape(-1, 1)
        if float(p.T @ p_prev) < 0.0:
            p = -p                              # 이전 스텝과 연속인 쪽
    return p
