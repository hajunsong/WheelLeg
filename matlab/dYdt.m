function [Yp, out] = dYdt(t, Y, prm)
% DYDT  상태벡터 Y 를 받아 미분 Yp 를 돌려준다. (RK4 의 각 stage 에서 호출)
%
%   Y  = [ r0(3) ; p0(4) ; q1 ;  dr0(3) ;  w0(3) ; dq1  ]   (15x1)
%   Yp = [ dr0(3); dp0(4); dq1;  ddr0(3);  dw0(3); ddq1 ]   (15x1)
%
%   out : 후처리/검증용 구조체 (CM 위치, 가속도, 구속 후 M/Q 등)

    % ---------------- Y2qdq ----------------
    r0  = Y(1:3, 1);
    p0  = Y(4:7, 1);
    q1  = Y(8, 1);
    dr0 = Y(9:11, 1);
    w0  = Y(12:14, 1);
    dq1 = Y(15, 1);

    % RK4 중간 stage 에서 노름이 흐르므로 매번 투영해준다
    p0 = p0/norm(p0);

    % ---------------- base body ----------------
    e0 = p0(2:4, 1);
    E0 = [-e0,  tilde(e0) + p0(1)*eye(3)];
    G0 = [-e0, -tilde(e0) + p0(1)*eye(3)];
    A0 = E0*G0';

    J0c = A0*prm.C00*prm.J0p*(A0*prm.C00)';

    rho0 = A0*prm.rho0p;
    r0c  = r0 + rho0;

    w0t = tilde(w0);
    r0t = tilde(r0);

    dr0c = dr0 + w0t*rho0;

    dr0t  = tilde(dr0);
    dr0ct = tilde(dr0c);
    r0ct  = tilde(r0c);

    Y0h = [dr0 + r0t*w0; w0];

    % 외력 : RecurDyn TRANSLATIONAL_FORCE, FY = step5(time, 0, F, 1, -F), RM = global
    %        작용점 base.Marker4 가 base.CM 과 같은 위치라 순수 CM 힘으로 들어간다
    Fy = step5(t, 0, prm.F_ex, 1, -prm.F_ex);

    f0c = [0; Fy; prm.m0*prm.g];
    t0c = [0; 0; 0];

    M0h = [prm.m0*eye(3), -prm.m0*r0ct;
           prm.m0*r0ct,    J0c - prm.m0*r0ct*r0ct];
    Q0h = [f0c + prm.m0*dr0ct*w0;
           t0c + r0ct*f0c + prm.m0*r0ct*dr0ct*w0 - w0t*J0c*w0];

    % ---------------- link body ----------------
    A01pp = [cos(q1), -sin(q1), 0;
             sin(q1),  cos(q1), 0;
             0,        0,       1];
    A1  = A0*prm.C01*A01pp;
    s01 = A1*prm.s01p;
    r1  = r0 + s01;

    J1c = A1*prm.C11*prm.J1p*(A1*prm.C11)';

    rho1 = A1*prm.rho1p;
    r1c  = r1 + rho1;
    r1ct = tilde(r1c);

    H1  = A0*prm.C01*[0;0;1];
    w1  = w0 + H1*dq1;
    w1t = tilde(w1);
    dr1 = dr0 + w0t*s01;
    r1t = tilde(r1);

    B1 = [r1t*H1; H1];

    dr1t  = tilde(dr1);
    dr1c  = dr1 + w1t*rho1;
    dr1ct = tilde(dr1c);

    dH1 = w0t*H1;
    D1  = [dr1t*H1 + r1t*dH1; dH1]*dq1;
    Y1h = Y0h + B1*dq1;

    f1c = [0; 0; prm.m1*prm.g];
    t1c = [0; 0; 0];

    M1h = [prm.m1*eye(3), -prm.m1*r1ct;
           prm.m1*r1ct,    J1c - prm.m1*r1ct*r1ct];
    Q1h = [f1c + prm.m1*dr1ct*w1;
           t1c + r1ct*f1c + prm.m1*r1ct*dr1ct*w1 - w1t*J1c*w1];

    % ---------------- mass / force ----------------
    K1 = M1h;
    K0 = K1 + M0h;

    L1 = Q1h;
    L0 = L1 + Q0h - K1*D1;

    % ---------------- EQM ----------------
    M = [K0,      K1*B1;
         B1'*K1,  B1'*K1*B1];
    Q = [L0;
         B1'*(L1 - K1*D1)];

    % 구속 : prm.free 이외의 자유도는 가속도 0 으로 잠근다
    lock = setdiff(1:7, prm.free);
    M(lock, :) = 0;
    M(sub2ind([7 7], lock, lock)) = 1;
    Q(lock) = 0;

    ddq  = M\Q;
    dY0h = ddq(1:6, 1);
    ddq1 = ddq(7, 1);

    % ---------------- base body acceleration ----------------
    dp0 = 0.5*E0'*w0;

    T0 = [eye(3),   -r0t;
          zeros(3),  eye(3)];
    R0 = [dr0t*w0; zeros(3,1)];
    dY0b = T0*dY0h - R0;

    ddr0 = dY0b(1:3, 1);
    dw0  = dY0b(4:6, 1);

    dw0t  = tilde(dw0);
    ddr0c = ddr0 + dw0t*rho0 + w0t*w0t*rho0;

    % ---------------- link body acceleration ----------------
    dY1h = dY0h + B1*ddq1 + D1;

    T1  = [eye(3),   -r1t;
           zeros(3),  eye(3)];
    dT1 = [zeros(3), -dr1t;
           zeros(3),  zeros(3)];
    dY1b = dT1*Y1h + T1*dY1h;

    ddr1 = dY1b(1:3, 1);
    dw1  = dY1b(4:6, 1);
    dw1t = tilde(dw1);

    ddr1c = ddr1 + dw1t*rho1 + w1t*w1t*rho1;

    % ---------------- dqddq2Yp ----------------
    Yp = [dr0; dp0; dq1; ddr0; dw0; ddq1];

    if nargout > 1
        out = struct('A0',A0, 'A1',A1, 'r0c',r0c, 'r1c',r1c, ...
                     'dr0c',dr0c, 'dr1c',dr1c, 'ddr0c',ddr0c, 'ddr1c',ddr1c, ...
                     'w1',w1, 'dw0',dw0, 'dw1',dw1, 'ddq1',ddq1, ...
                     'Fy',Fy, 'M',M, 'Q',Q);
    end
end
