function xdot = system_dynamics(t, x, u, current_field, wheelbase, Q, H, R)
% Compute time derivative of the full vehicle and covariance state.
%
% Inputs:
%   t             - time [s]
%   x             - full state (39 x n) [pE; pN; theta; L_vec]
%   u             - control input (2 x n) [v; delta] [m/s; rad]
%   current_field - function handle mapping (t, x, y) to (2 x n) current velocity [m/s]
%   wheelbase     - effective vehicle length governing turn rate in bicycle kinematic model [m]
%   Q             - INS process noise covariance matrix (8 x 8)
%   H             - measurement Jacobian (m x 8)
%   R             - measurement noise covariance matrix (m x m)
%
% Outputs:
%   xdot - time derivative of full state (39 x n)

xdot = [nav.uuv_dynamics(t, x(1:3,:), u, current_field, wheelbase);
        nav.cov_dynamics(-x(3,:), x(4:end,:), Q, H, R)];
end
