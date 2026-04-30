function [Lpdot_vec] = pos_cov_dynamics(pos, Lp_vec, Lvss_vec, anom_meas)
% Compute time derivative of the Cholesky factor of the position error covariance
% using the quasi-static approximation with corections based on geophysical anomalies.
%
% Inputs:
%   pos       - vehicle position (2 x n) [pE; pN] [m]
%   Lp_vec    - lower Cholesky factor of position error covariance, vectorized (3 x n) [m]
%   Lvss_vec  - lower Cholesky factor of steady-state velocity error covariance, vectorized (3 x 1) [m/s]
%   anom_meas - function handle mapping (x, y) to H_anom (2 x 2) and R_anom (2 x 2)
%
% Outputs:
%   Lpdot_vec - time derivative of lower Cholesky factor of position error covariance, vectorized (3 x n) [m/s]

n         = size(Lp_vec, 2);
Lpdot_vec = zeros(3, n);

Lvss       = utils.lvec_to_L(Lvss_vec);
Pdot_drift = 2 * (Lvss * Lvss');

for k = 1:n
    [H_anom, R_anom] = anom_meas(pos(1,k), pos(2,k));

    L    = utils.lvec_to_L(Lp_vec(:,k));
    Pp   = L * L';
    Pdot = Pdot_drift - Pp * H_anom' * (R_anom \ (H_anom * Pp));

    M               = L \ Pdot / L';
    S               = tril(M, -1) + 0.5 * diag(diag(M));
    Lpdot           = L * S;
    Lpdot_vec(:, k) = utils.L_to_lvec(Lpdot);
end
end
