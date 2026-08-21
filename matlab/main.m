clc; clear; close all;

%% ======================= parameter =======================
% 단위계 : SI (m, kg, s, N). 관성모멘트 단위는 kg*m^2 이다.

% ---- base body (기준계 = base.Ai) ----
prm.rho0p = [0;0;0];
prm.C00   = ang2mat(0, pi/2, 0);

prm.m0 = 111.015764646288;
Ixx = 9.34382685772924;  Ixy = 0;
Iyy = 0.27753941161572;  Iyz = 0;
Izz = 9.43633999493448;  Izx = 0;
prm.J0p = [Ixx, Ixy, Izx;
           Ixy, Iyy, Iyz;
           Izx, Iyz, Izz];

% ---- link body (기준계 = body.Ai) ----
prm.s01p  = [0;0;0];
prm.C01   = ang2mat(pi/2, pi/2, pi/2);

prm.rho1p = [6.22253616772498e-6; 0.359132273856264; 0];  % body.CM QP - body.Ai QP [m] (Ai 프레임)
prm.C11   = ang2mat(pi, pi/2, pi/2);

prm.m1 = 18.6342592049469;
Ixx = 0.604104103452763;  Ixy =  1.26126031622254e-7;
Iyy = 0.56572829220829;   Iyz = -3.57293476271758e-5;
Izz = 0.0403172972174389; Izx = -3.51625272013543e-19;
prm.J1p = [Ixx, Ixy, Izx;
           Ixy, Iyy, Iyz;
           Izx, Iyz, Izz];

% ---- system ----
prm.g    = -9.80665;        % m/s^2
prm.F_ex = 10;              % cart_pole_16 : FY = step5(time, 0, 10, 1, -10) [N]
prm.free = [2, 7];          % 살릴 자유도 : base 전역 Y 병진 + 회전 조인트

h   = 0.001;
t_e = 2;

%% ======================= initial condition =======================
r0  = [0;0;-0.05];                      % base.Ai 의 전역 위치 [m]
p0  = mat2ep(ang2mat(0, -pi/2, 0));     % base.Ai 의 전역 자세
q1  = 0;
dr0 = [0;0;0];
w0  = [0;0;0];
dq1 = 3;                                % 주의: RecurDyn cart_pole_16 은 초기속도 0 이라
                                        %       아래 RecurDyn 비교표는 이 값에서 어긋난다

Y = [r0; p0; q1; dr0; w0; dq1];

%% ======================= RK4 =======================
n  = round(t_e/h);
T  = (0:n)*h;               % 부동소수 누적 오차 없이 정확한 격자
YY = zeros(15, n+1);
AA = zeros(15, n+1);        % Yp 이력 (가속도 비교용)

YY(:,1) = Y;

for k = 1:n
    t  = T(k);
    k1 = dYdt(t,       Y,          prm);
    k2 = dYdt(t + h/2, Y + h/2*k1, prm);
    k3 = dYdt(t + h/2, Y + h/2*k2, prm);
    k4 = dYdt(t + h,   Y + h*k3,   prm);

    AA(:,k) = k1;
    Y = Y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

    Y(4:7) = Y(4:7)/norm(Y(4:7));       % 오일러 파라미터 정규화

    YY(:,k+1) = Y;
end
AA(:,end) = dYdt(T(end), YY(:,end), prm);

%% ======================= post processing =======================
% cart_pole_13 의 cart_px/py/pz = dx/dy/dz(base.Ai, Ground.origin, Ground.origin)
% Ground.origin 은 QP=(0,0,0), REULER=(0,0,0) 즉 전역계라 변환이 필요없다.
cart_p = YY(1:3,  :);       % [cart_px; cart_py; cart_pz]  [m]
cart_v = YY(9:11, :);       % [cart_vx; cart_vy; cart_vz]  [m/s]
cart_a = AA(9:11, :);       % [cart_accx; ...]             [m/s^2]

pend_q   = YY(8,  :);       % pendulum_q   = az(body.Ai, base.Cij)   [rad]
pend_qd  = YY(15, :);       % pendulum_qd                            [rad/s]
pend_qdd = AA(15, :);       % pendulum_qdd                           [rad/s^2]

fprintf('t = %.3f ~ %.3f s, %d steps (h = %g)\n', T(1), T(end), n, h);
fprintf('quaternion norm drift : %.3e\n', max(abs(vecnorm(YY(4:7,:)) - 1)));

%% ======================= RecurDyn 비교 =======================
% rec_data.csv : 헤더 없음, 8 열 (RecurDyn 원본 병진 채널은 mm 단위)
%   [ index, time, cart_py, cart_vy, cart_accy, pendulum_q, pendulum_qd, pendulum_qdd ]
csv = fullfile(fileparts(mfilename('fullpath')), '..', 'recurdyn', '01_cart_pole', 'rec_data.csv');

lab  = {'cart\_py [m]',     'cart\_vy [m/s]',      'cart\_accy [m/s^2]', ...
        'pendulum\_q [rad]','pendulum\_qd [rad/s]','pendulum\_qdd [rad/s^2]'};
mine = [cart_p(2,:); cart_v(2,:); cart_a(2,:); pend_q; pend_qd; pend_qdd];

ref = [];
if isfile(csv)
    R     = readmatrix(csv);
    ref.t = R(:,2).';
    ref.y = R(:,3:8).';
    ref.y(1:3,:) = 1e-3*ref.y(1:3,:);              % mm 계열 -> SI(m 계열)
    fprintf('\nRecurDyn : rec_data.csv (%d points, t = %.3f ~ %.3f)\n', ...
            numel(ref.t), ref.t(1), ref.t(end));
    fprintf('%-24s %12s %12s %10s\n', 'channel', 'max|err|', 'RMS', 'rel.RMS');
    err = zeros(6, numel(ref.t));
    for i = 1:6
        ip       = interp1(T, mine(i,:), ref.t, 'linear');
        err(i,:) = ip - ref.y(i,:);
        sc       = max(abs(ref.y(i,:)));
        fprintf('%-24s %12.4e %12.4e %10.2e\n', erase(lab{i},'\'), ...
                max(abs(err(i,:))), rms(err(i,:)), rms(err(i,:))/sc);
    end
    % RecurDyn 의 acc 채널은 알고리즘 감쇠(NDAMPING)로 오염되어 있다.
    % 자기 속도를 미분한 값과 비교해야 실제 모델 일치도가 보인다.
    dt  = ref.t(2) - ref.t(1);
    anm = {'cart_accy', 'pendulum_qdd'};
    for j = 1:2
        cc = [1 4]; c = cc(j);
        a_num = gradient(ref.y(c+1,:), dt);              % RecurDyn 의 d(vel)/dt
        kk    = true(size(ref.t));                       % step5 라 불연속 없음
        am    = interp1(T, mine(c+2,:), ref.t, 'linear');
        sc    = max(abs(ref.y(c+2,:)));
        fprintf('  %-14s vs 보고 acc %8.2e | vs d(vel)/dt %8.2e  (rel.RMS)\n', ...
                anm{j}, rms(am(kk) - ref.y(c+2,kk))/sc, rms(am(kk) - a_num(kk))/sc);
    end
else
    fprintf('\n[비교 생략] %s 없음\n', csv);
end

%% ======================= plot =======================
figure('Name','MATLAB vs RecurDyn','Position',[60 60 1150 750]);
ord = [1 3 5 2 4 6];                        % 좌열 cart, 우열 pendulum
for i = 1:6
    subplot(3,2,ord(i));
    plot(T, mine(i,:), 'b', 'LineWidth', 1.4); hold on;
    if ~isempty(ref), plot(ref.t, ref.y(i,:), 'r--', 'LineWidth', 1.2); end
    grid on; xlabel('time [s]'); ylabel(lab{i});
    if i == 1 && ~isempty(ref), legend('MATLAB','RecurDyn','Location','best'); end
end

if ~isempty(ref)
    figure('Name','error (MATLAB - RecurDyn)','Position',[100 100 1150 750]);
    for i = 1:6
        subplot(3,2,ord(i));
        plot(ref.t, err(i,:), 'k', 'LineWidth', 1.1);
        grid on; xlabel('time [s]'); ylabel(['\Delta ' lab{i}]);
    end
end
