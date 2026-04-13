function [Lpdot_vec] = pos_cov_dynamics(Lp_vec, Lvss_vec)
% Compute time derivative of the Cholesky factor of the position error covariance
% using the quasi-static approximation: Pdot_pos ≈ P_vel_ss = Lvss * Lvss'.
%
% Inputs:
%   Lp_vec   - lower Cholesky factor of position error covariance, vectorized (3 x n) [m]
%   Lvss_vec - lower Cholesky factor of steady-state velocity error covariance, vectorized (3 x 1) [m/s]
%
% Outputs:
%   Lpdot_vec - time derivative of lower Cholesky factor of position error covariance, vectorized (36 x n) [m/s]

n         = size(Lp_vec, 2);
Lpdot_vec = zeros(3, n);

Lvss = utils.lvec_to_L(Lvss_vec);
Pdot = 2 * (Lvss * Lvss');

for k = 1:n
    L               = utils.lvec_to_L(Lp_vec(:,k));
    M               = L \ Pdot / L';
    S               = tril(M, -1) + 0.5 * diag(diag(M));
    Lpdot           = L * S;
    Lpdot_vec(:, k) = utils.L_to_lvec(Lpdot);
end
end
