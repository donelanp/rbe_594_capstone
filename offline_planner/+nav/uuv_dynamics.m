function [xdot] = uuv_dynamics(t, x, u, current_field, wheelbase)
% Compute time derivative of vehicle state under bicycle kinematic model with ocean current.
%
% Inputs:
%   t             - time [s]
%   x             - vehicle state (3 x n) [pE; pN; theta] [m; m; rad]
%   u             - control input (2 x n) [v; delta] [m/s; rad]
%   current_field - function handle mapping (t, x, y) to (2 x n) current velocity [m/s]
%   wheelbase     - effective vehicle length governing turn rate in bicycle kinematic model [m]
%
% Outputs:
%   xdot - time derivative of vehicle state (3 x n) [pE_dot; pN_dot; theta_dot] [m/s; m/s; rad/s]

theta = x(3,:);
v     = u(1,:);
delta = u(2,:);
c     = current_field(t, x(1,:), x(2,:));

xdot = [v .* sin(theta) + c(1,:);
        v .* cos(theta) + c(2,:);
        (v / wheelbase) .* tan(delta)];
end
