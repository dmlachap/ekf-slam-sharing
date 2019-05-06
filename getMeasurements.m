function z = getMeasurements(x, Q, c, landmarks, phirange)
% landmarks: Nx2 matrix of x and y values of landmarks
% phirange: range of valid bearings for visible landmarks

N = length(c);

% indices in z of landmark i: (3*i-2) : (3*i)
z = zeros(3*N, 1);
for i = 1:N
    phi = atan2(landmarks(i,2) - x(2), landmarks(i,1)- x(1)) - x(3);
%     phi = wrapTo2Pi(phi); % bearing angle from 0 (directly ahead)
    if wrapTo2Pi(phi) >= wrapTo2Pi(phirange(1)) || wrapTo2Pi(phi) <= phirange(2)
        z(3*i-2:3*i) = [norm(landmarks(i,:) - x(1:2).') + normrnd(0, Q(1,1));
                        wrapTo2Pi(phi + normrnd(0, Q(2,2)));
                        i];
    end
end


end