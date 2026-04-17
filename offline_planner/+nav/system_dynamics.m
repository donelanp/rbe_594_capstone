function [xdot] = system_dynamics(t, x, u, current_field, wheelbase, Lvss_vec)
% Compute time derivative of the full vehicle and covariance state.
%
% Inputs:
%   t             - time [s]
%   x             - full state (6 x n) [pE; pN; theta; Lp_vec]
%   u             - control input (2 x n) [v; delta] [m/s; rad]
%   current_field - function handle mapping (t, x, y) to (2 x n) current velocity [m/s]
%   wheelbase     - effective vehicle length governing turn rate in bicycle kinematic model [m]
%   Lvss_vec      - lower Cholesky factor of steady-state velocity error covariance, vectorized (3 x 1) [m/s]
%
% Outputs:
%   xdot - time derivative of full state (6 x n)

xdot = [nav.uuv_dynamics(t, x(1:3,:), u, current_field, wheelbase);
        nav.pos_cov_dynamics(x(4:6,:), Lvss_vec)];
end
