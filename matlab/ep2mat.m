function A = ep2mat(p)
% EP2MAT  오일러 파라미터 [e0;e1;e2;e3] -> 회전행렬.  mat2ep 의 역변환.

    p  = p(:)/norm(p);
    e0 = p(1);
    e  = p(2:4);

    E = [-e,  tilde(e) + e0*eye(3)];
    G = [-e, -tilde(e) + e0*eye(3)];
    A = E*G';
end
