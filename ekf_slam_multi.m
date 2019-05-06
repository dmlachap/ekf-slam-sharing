% ekf_slam_multi
% close all;
% clear;

%% Environment Setup
K = 10; % number of robots

dt = 0.1; % time step
T = 20; % end time
tvec = 0:dt:T;

rng(3);

bound = 5; % environment boundary
axlimits = [-bound bound -bound bound];
N = 12;
landmarks = 2*bound*rand(N, 2) - bound; % landmark coordinates

c = 1:N; % landmark signatures / correspondence variables

phirange = [-pi/4 pi/4]; % maximum sense-able landmark bearing

limit = N;

%% Dynamics
syms x y theta v w
xsymvec = [x y theta].';
usymvec = [v w].';

xnew = x -(v/w)*sin(theta) + (v/w)*sin(theta + w*dt);
ynew = y + (v/w)*cos(theta) - (v/w)*cos(theta + w*dt);
thetanew = theta + w*dt;

fsym = [xnew;
        ynew;
        thetanew];
Fsym = jacobian(fsym, xsymvec);

f = matlabFunction(fsym, 'Vars', {xsymvec, usymvec});
F = matlabFunction(Fsym, 'Vars', {xsymvec, usymvec});

%% Initialize Robots
nt = length(tvec);
nx = length(fsym);
robots = struct;
for k = 1:K
    theta0 = 2*pi*rand(1);
    Rot = [cos(theta0) -sin(theta0);
         sin(theta0) cos(theta0)];
     
    robots(k).x0G = [2*bound*rand(2, 1) - bound; theta0]; % initial state (in global frame)
    robots(k).Rot = Rot; % rotation matrix
    % note that each robot *thinks* that it is at origin, pointed along
    % positive x-axis initially
    
    % true state
    robots(k).xtrue = zeros(nx, nt);
    
    % measurements
    robots(k).z = zeros(3*N, nt);
    
    % state estimate
    robots(k).xhat = zeros(nx + 3*N, nt);
    
    % error covariance estimate
    robots(k).P = cell(1, nt);
    P0 = 100*eye(3 + 3*N); P0(1:3, 1:3) = 0;
    robots(k).P{1} = P0;
    
    % control inputs
    robots(k).u = [abs(normrnd(0.1, 0.5))*ones(1, nt);
                   sign(rand(1)-0.5)*normrnd(0.5, 0.1)*ones(1, nt)];
    
    robots(k).seen = false(1, N); % keeps track of seen landmarks
    
    robots(k).R = 0.0001*eye(3); % process noise covariance
    
    robots(k).Q = 0.001*eye(3); % measurement noise covariance
    
end

% map-sharing struct
transmit = struct;
transmit.maps = cell(K, nt); % vector of maps to transmit
transmit.P = cell(K, nt);

% performance measures
logdetP = zeros(K, nt);
state_squared_error = zeros(K, nt);
map_squared_error = zeros(K, nt);
nees = zeros(K, nt);

%% Initialize video
v = VideoWriter(sprintf('EKF SLAM K=%d L=%d.avi', K, limit), 'Uncompressed AVI');
v.FrameRate = round(1/dt);
open(v);

%% Run Filter
figure('Position', [0 0 1024 1024])
ax = axes;
cmap = lines;
for ti = 2:nt % loop over time steps
    tic
    cla;
    axis(ax, axlimits)
    axis tight equal
    grid on
    hold on
    t = tvec(ti);
    
    for k = 1:K % loop over robots
        color = cmap(k, :);
        
        gMk = [robots(k).Rot robots(k).x0G(1:2);
               0 0 1];
        
        % true dynamics (local frame)
        robots(k).xtrue(:, ti) = f(robots(k).xtrue(:, ti-1), robots(k).u(:, ti-1)) + mvnrnd(zeros(1, nx), robots(k).R).';
        
        % Transform landmarks into local frame
        landmarks_k = (gMk\([landmarks ones(N, 1)].')).';
        landmarks_k(:, 3) = 1:N;
        landmarks_k_xy = landmarks_k(:, 1:2);
        landmarks_k_flat = landmarks_k.'; landmarks_k_flat = landmarks_k_flat(:);
        robots(k).z(:, ti) = getMeasurements(robots(k).xtrue(:, ti), robots(k).Q, c, landmarks_k_xy, phirange);
        
        % ekf update
        [robots, transmit] = ekf_slam_known_sharing(ti, k, robots, f, F, c, transmit, limit);
        
        % draw environment
        drawEnvironment(ax, robots(k).xhat(:, ti), robots(k).P{ti}, robots(k).xtrue(:, ti), landmarks, color, gMk);
        
        % performance
        logdetP(k, ti) = log(det(robots(k).P{ti}));
        delta = [robots(k).xtrue(:, ti); landmarks_k_flat] - robots(k).xhat(:, ti);
        state_squared_error(k, ti) = delta(1:3).'*delta(1:3);
        map_delta = delta(4:end); map_delta = map_delta(robots(k).xhat(4:end,ti)~=0);
        map_squared_error(k, ti) = map_delta.'*map_delta/length(map_delta);
        nees(k, ti) = delta(1:3).'*inv(robots(k).P{ti}(1:3, 1:3))*delta(1:3);
    end
    
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    hold off;
    pause(dt-toc);   
end

close(v);

%% Results
figure()
hold on;
for k = 1:K
    plot(tvec, logdetP(k,:))
end
title(sprintf('Estimated Log-Determinant Covariance (Full State), K=%d, L=%d', K, limit))
xlabel('Time (s)')
ylabel('logdet(P)')
hold off;

figure()
hold on;
for k = 1:K
    plot(tvec, state_squared_error(k,:))
end
title(sprintf('Squared Error (Robot State), K=%d, L=%d', K, limit))
xlabel('Time (s)')
ylabel('x_R^Tx_R')
hold off;

figure()
hold on;
for k = 1:K
    plot(tvec, map_squared_error(k,:))
end
title(sprintf('Squared Error (Map State), K=%d, L=%d', K, limit))
xlabel('Time (s)')
ylabel('x_M^Tx_M')
hold off;

figure()
hold on;
alpha = 0.05;
lowbound = chi2inv(alpha/2, nx);
highbound = chi2inv(1 - alpha/2, nx);
for k = 1:K
    plot(tvec, nees(k, :))
end
plot([tvec(1) tvec(end)], [lowbound, lowbound], 'r--', ...
    [tvec(1) tvec(end)], [highbound, highbound], 'r--')
title(sprintf('Robot NEES, K=%d, L=%d', K, limit))
xlabel('Time (s)')
ylabel('x_R^TP_R^{-1}x_R')
hold off;