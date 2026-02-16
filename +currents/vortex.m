function [v_c] = vortex(x, y, x_center, y_center, intensity)
% Compute vortex current centered at a point.
%
% Inputs:
%   x         - x positions [m]
%   y         - y positions [m]
%   x_center  - vortex center x coordinate [m]
%   y_center  - vortex center y coordinate [m]
%   intensity - vortex strength [m^2/s]
%
% Outputs:
%   v_c - current velocity (2 x n) [m/s]

dx      = x - x_center;
dy      = y - y_center;
r2      = dx.^2 + dy.^2 + 1e-6;
v_theta = intensity ./ sqrt(r2);
v_c     = [-v_theta .* dy ./ sqrt(r2);
            v_theta .* dx ./ sqrt(r2)];
end
