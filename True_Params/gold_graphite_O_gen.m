clear; clc; close all; hold on;

rng(1357354148)

v = noisy_y(deg2rad(20), 10);

[X,Y,Z] = sphere(50);   % 50 controls resolution
surf(X,Y,Z,FaceColor='b', FaceAlpha=0.2, EdgeColor='none')
plot3([0,1],[0,0],[0,0], 'r')
plot3([0,0],[0,1],[0,0], 'g')
plot3([0,0],[0,0],[0,1], 'b')
quiver3(zeros(size(v(:,1))),zeros(size(v(:,1))),zeros(size(v(:,1))),v(:,1), v(:,2), v(:,3), 1, Color='k')
axis equal

save("gold_graphite_O.mat", "v")

function v = noisy_y(sigma, N)    
    % Sample polar angle (angle from +y axis)
    theta = normrnd(0, sigma, N, 1);   % small-angle Gaussian
    
    % Sample azimuth uniformly
    phi = pi*rand(N,1);
    
    % Convert to Cartesian (y is the "pole")
    x = sin(theta) .* cos(phi);
    y = cos(theta);
    z = sin(theta) .* sin(phi);
    
    % Combine into unit vectors
    v = [x(:) y(:) z(:)];
end