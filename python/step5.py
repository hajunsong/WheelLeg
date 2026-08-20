"""step5.m 1:1 변환.  ADAMS/RecurDyn 의 5차 다항 step 함수."""


def step5(x, x0, h0, x1, h1):
    """
    x <= x0 : h0
    x >= x1 : h1
    그 사이 : h0 + (h1-h0)*a^3*(10 - 15a + 6a^2),   a = (x-x0)/(x1-x0)

    양 끝에서 1계, 2계 도함수가 모두 0 이라 C2 연속이다.
    (3차 STEP 은 C1 까지만 연속이라 가속도가 꺾인다)
    """
    if x <= x0:
        y = h0
    elif x >= x1:
        y = h1
    else:
        a = (x - x0) / (x1 - x0)
        y = h0 + (h1 - h0) * a**3 * (10.0 - 15.0*a + 6.0*a**2)
    return y
