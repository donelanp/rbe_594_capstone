function [xdot] = system_dynamics(t, x, u, current_field, wheelbase, Lvss_vec, anom_meas)
% Compute time derivative of the full vehicle and covariance state.
%
% Inputs:
%   t             - time [s]
%   x             - full state (6 x n) [pE; pN; theta; Lp_vec]
%   u             - control input (2 x n) [v; delta] [m/s; rad]
%   current_field - function handle mapping (t, x, y) to (2 x n) current velocity [m/s]
%   wheelbase     - effective vehicle length governing turn rate in bicycle kinematic model [m]
%   Lvss_vec      - lower Cholesky factor of steady-state velocity error covariance, vectorized (3 x 1) [m/s]
%   anom_meas     - function handle mapping (x, y) to H_anom (2 x 2) and R_anom (2 x 2)
%
% Outputs:
%   xdot - time derivative of full state (6 x n)

xdot = [nav.uuv_dynamics(t, x(1:3,:), u, current_field, wheelbase);
        nav.pos_cov_dynamics(x(1:2,:), x(4:6,:), Lvss_vec, anom_meas)];
end
