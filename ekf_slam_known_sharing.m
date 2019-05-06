function [robots, transmit] = ekf_slam_known_sharing(ti, robot_idx, robots, f, F, c, transmit, limit)

N = length(c);
Ix = [eye(3) zeros(3, 3*N)];

% unpack robot structure
x = robots(robot_idx).xhat(:, ti-1);
P = robots(robot_idx).P{ti-1};
u = robots(robot_idx).u(:, ti);
z = robots(robot_idx).z(:, ti);
R = robots(robot_idx).R;
Q = robots(robot_idx).Q;
seen = robots(robot_idx).seen;

xrobot = x(1:3);
xmap = [zeros(3, 1); x(4:end)];
xbar = Ix.'*f(xrobot, u) + xmap;

Gt = eye(length(xbar)) + Ix.'*(F(xrobot, u)-eye(3))*Ix;

Pbar = Gt*P*Gt.' + Ix.'*R*Ix;

Hifull = [];
zifull = [];
zhatifull = [];
Qfull = [];
seenthistime = 0; % number of landmarks seen
for i = 1:N % loop over landmarks
    if z(3*i)~=0 % valid landmark seen
        seenthistime = seenthistime + 1;
        zi = z(3*i-2:3*i);
        if ~seen(i) % initialize landmark position
            xbar(3 + 3*i-2 : 3 + 3*i) = [xbar(1) + zi(1)*cos(zi(2) + xbar(3));
                                         xbar(2) + zi(1)*sin(zi(2) + xbar(3));
                                         zi(3)];
            seen(i) = true;
        end
        
        % predicted measurements
        delta = [xbar(3 + 3*i-2) - xbar(1);
                 xbar(3 + 3*i-1) - xbar(2)];
        q = delta.'*delta;
        zi(2) = wrapTo2Pi(zi(2) - pi);
        phihat = wrapTo2Pi(atan2(delta(2), delta(1)) - xbar(3) - pi);
        zhati = [sqrt(q);
                 phihat;
                 xbar(3 + 3*i)];
             
        Ixi = [eye(3) zeros(3, 3*N);
               zeros(3, 3*i) eye(3) zeros(3, 3*N - 3*i)];
        
        Hifull = [Hifull; 
            (1/q)*[-sqrt(q)*delta(1) -sqrt(q)*delta(2) 0 sqrt(q)*delta(1) sqrt(q)*delta(2) 0;
                     delta(2) -delta(1) -q -delta(2) delta(1) 0;
                     0 0 0 0 0 q]*Ixi];
        zifull = [zifull; 
                  zi];
        zhatifull = [zhatifull;
                     zhati];
        Qfull = blkdiag(Qfull, Q);
    end
    
end

% rotation matrix, global<-robot local
gMe = [robots(robot_idx).Rot robots(robot_idx).x0G(1:2);
       0 0 1];

shared_maps = transmit.maps(:, ti-1);
shared_P = transmit.P(:, ti-1);
for k = 1:length(shared_maps) % loop over robots (to analyze maps)
    if k ~= robot_idx
        kmap = shared_maps{k};
        kP = shared_P{k};
        for j = 1:length(kmap)/3 % loop over landmarks in map s
            k_localmap = kmap(3*j - 2 : 3*j);

            gMk = [robots(k).Rot robots(k).x0G(1:2);
                    0 0 1];
            
            e_localmap = gMe\gMk*[k_localmap(1:2); 1];

            zj = [e_localmap(1:2); k_localmap(3)];

            i = zj(3);
            if ~seen(i)
               xbar(3 + 3*i-2 : 3+ 3*i) = zj;
               seen(i) = true; 
            else
                Hifull = [Hifull; 
                          zeros(3, 3*i) eye(3) zeros(3, 3*N - 3*i)];
                zifull = [zifull; 
                          zj];
                zhatifull = [zhatifull;
                             xbar(3 + 3*i - 2 : 3 + 3*i)];
                
                b = 0.1;
                Qfull = blkdiag(Qfull, kP{j} + b*R);
            end
        end
    end
end

if seenthistime > 0
    Ki = Pbar*Hifull.'/(Hifull*Pbar*Hifull.' + Qfull);
    xbar = xbar + Ki*(zifull - zhatifull);
    Pbar = (eye(size(Pbar)) - Ki*Hifull)*Pbar;
end

% pack robot structure
robots(robot_idx).xhat(:, ti) = xbar;
robots(robot_idx).P{ti} = Pbar;
robots(robot_idx).seen = seen;

% determine map fragment to share
valid_map_pts = xbar(4:end); valid_map_pts = valid_map_pts(valid_map_pts ~=0);
valid_map_idxs = valid_map_pts(3:3:end);
shareorder = [];
for m = valid_map_idxs.'
    shareorder = [shareorder; -log(det(Pbar(3*m - 2 : 3*m, 3*m - 2 : 3*m))), m];
end
shareorder = sortrows(shareorder);
toshare = [];
shareP = {};
if ~isempty(shareorder)
    l = 1;
    while l < limit && l < length(valid_map_idxs)
        shareidx = shareorder(l, 2);
        toshare = [toshare; xbar(3 + 3*shareidx - 2 : 3 + 3*shareidx)];
        shareP{end+1} = Pbar(3 + 3*shareidx - 2 : 3 + 3*shareidx, 3 + 3*shareidx - 2 : 3 + 3*shareidx);
        l = l+1;
    end
end

% pack transmit structure
transmit.maps{robot_idx, ti} = toshare;
transmit.P{robot_idx, ti} = shareP;

end