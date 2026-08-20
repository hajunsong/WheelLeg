function p = mat2ep(A, p_prev)
% MAT2EP  회전행렬 -> 오일러 파라미터 (단위 쿼터니언, scalar first)
%
%   p = mat2ep(A)          A(3x3) -> p = [e0;e1;e2;e3],  e0 >= 0 으로 정규화
%   p = mat2ep(A, p_prev)  p_prev 와 부호가 연속이 되도록 선택 (적분 루프용)
%
%   규약 (main.m 의 E, G 와 동일):
%       E = [-e,  tilde(e) + e0*eye(3)]
%       G = [-e, -tilde(e) + e0*eye(3)]
%       A = E*G'
%
%   Shepperd 방법: 4*ei^2 후보 중 가장 큰 것부터 구해 0 나눗셈을 피한다.

    tr = A(1,1) + A(2,2) + A(3,3);

    % 4*e0^2, 4*e1^2, 4*e2^2, 4*e3^2
    c = [1 + tr;
         1 + 2*A(1,1) - tr;
         1 + 2*A(2,2) - tr;
         1 + 2*A(3,3) - tr];

    [cmax, k] = max(c);
    s = sqrt(cmax);     % 2*|e_k|
    d = 0.5/s;          % 1/(4*e_k)

    switch k
        case 1
            e0 = 0.5*s;
            e1 = (A(3,2) - A(2,3))*d;
            e2 = (A(1,3) - A(3,1))*d;
            e3 = (A(2,1) - A(1,2))*d;
        case 2
            e1 = 0.5*s;
            e0 = (A(3,2) - A(2,3))*d;
            e2 = (A(1,2) + A(2,1))*d;
            e3 = (A(1,3) + A(3,1))*d;
        case 3
            e2 = 0.5*s;
            e0 = (A(1,3) - A(3,1))*d;
            e1 = (A(1,2) + A(2,1))*d;
            e3 = (A(2,3) + A(3,2))*d;
        otherwise
            e3 = 0.5*s;
            e0 = (A(2,1) - A(1,2))*d;
            e1 = (A(1,3) + A(3,1))*d;
            e2 = (A(2,3) + A(3,2))*d;
    end

    p = [e0; e1; e2; e3];
    p = p/norm(p);

    % p 와 -p 는 같은 자세이므로 부호를 하나 골라야 한다
    if nargin < 2
        if p(1) < 0
            p = -p;                 % e0 >= 0  (theta in [0, 180deg])
        end
    else
        if p.'*p_prev(:) < 0
            p = -p;                 % 이전 스텝과 연속인 쪽
        end
    end
end
