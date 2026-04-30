function [Lvss_vec] = vel_cov_steady_state(Q, H, R)
% Compute the Cholesky factor of the steady-state velocity error covariance by
% solving the algebraic Riccati equation for the reduced observable state
% [dv_N, dv_E, dpsi, bax, bay, bgz] at yaw = 0.
%
% Inputs:
%   Q - INS process noise covariance (8 x 8)
%   H - measurement Jacobian (m x 8)
%   R - measurement noise covariance (m x m)
%
% Outputs:
%   Lvss_vec - lower Cholesky factor of steady-state velocity error covariance, vectorized (3 x 1) [m/s]

Q_r = Q(3:8, 3:8);
H_r = H(:, 3:8);

F_r       = zeros(6);
F_r(1, 4) = 1;
F_r(2, 5) = 1;
F_r(3, 6) = 1;

P_r      = icare(F_r', H_r', Q_r, R);
Lvss     = chol(P_r(1:2, 1:2), 'lower');
Lvss_vec = utils.L_to_lvec(Lvss);
end
