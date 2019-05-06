function drawEnvironment(ax, x, P, xtrue, landmarks, color, varargin)

% clf;
if nargin == 6 % x given in global frame
    drawRobot(ax, x(1:3), P(1:3, 1:3), color);
    drawRobot(ax, xtrue(1:3), P(1:3, 1:3), 'g');
    drawLandmarks(ax, x(4:end), P(4:end, 4:end), landmarks, color);
    shg;
elseif nargin == 7 % transformation global<-local given
    T = varargin{1};
    xG = T*[x(1:2); 1]; xG = [xG(1:2); acos(T(1,1)) + x(3)];
    PG = [T(1:2, 1:2)*P(1:2, 1:2), zeros(2, 1);
          zeros(1, 2),  P(3, 3)];
    
    xtrueG = T*[xtrue(1:2); 1]; xtrueG = [xtrueG(1:2); acos(T(1,1)) + xtrue(3)];
    
    drawRobot(ax, xG(1:3), PG(1:3, 1:3), color);
    drawRobot(ax, xtrueG(1:3), PG(1:3, 1:3), 'g');
    
    map = x(4:end); map(3:3:end) = 1;
    mapK = [map(1:3:end).';
            map(2:3:end).';
            map(3:3:end).'];
        
    mapK = mapK(:, mapK(1,:)~=0);
              
    mapG = T*mapK;
    mapG = mapG(:);
    
    drawLandmarks(ax, mapG, P(4:end, 4:end), landmarks, color);
    shg;
else
    error('Invalid number of arguments');
end


end

function drawRobot(ax, x, P, color)
r = 0.1; % robot radius

% hold on

th = 0:pi/50:2*pi;
xunit = r*cos(th) + x(1);
yunit = r*sin(th) + x(2);
% hold on
plot(ax, xunit, yunit, 'Color', color)
plot(ax, [x(1) x(1) + r*cos(x(3))], [x(2) x(2) + r*sin(x(3))], 'Color', color);
% hold off

% hold off

end

function drawLandmarks(ax, x, P, landmarks, color)

x = x(abs(x)>=1e-4);

% hold on
plot(ax, landmarks(:,1), landmarks(:,2), 'g*');
plot(ax, x(1:3:end), x(2:3:end), 'Color', color, 'Marker', '*', 'LineStyle', 'none');

% hold off

end