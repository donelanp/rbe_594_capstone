function [Ldot_vec] = cov_dynamics(psi, pos, L_vec, Q, H, R, anom_meas)
% Compute time derivative of the Cholesky factor of the INS error covariance.
%
% Inputs:
%   psi       - vehicle yaw (1 x n) [rad]
%   pos       - vehicle position (2 x n) [pE; pN] [m]
%   L_vec     - lower Cholesky factor of error covariance, vectorized (36 x n)
%   Q         - process noise covariance matrix (8 x 8)
%   H         - non-geophysical anomaly measurement Jacobian (m x 8)
%   R         - non-geophysical anomaly measurement noise covariance matrix (m x m)
%   anom_meas - function handle mapping (x, y) to H_anom (2 x 2) and R_anom (2 x 2)
%
% Outputs:
%   Ldot_vec - time derivative of lower Cholesky factor of error covariance, vectorized (36 x n)

N        = size(L_vec, 2);
Ldot_vec = zeros(36, N);

for k = 1:N
    L = utils.lvec_to_L(L_vec(:, k));
    P = L * L';

    cp = cos(psi(k));
    sp = sin(psi(k));

    F       = zeros(8);
    F(1, 3) = 1;
    F(2, 4) = 1;
    F(3, 6) = cp;
    F(3, 7) = sp;
    F(4, 6) = -sp;
    F(4, 7) = cp;
    F(5, 8) = 1;

    % augment observation model with geophysical anomaly measurements
    [H_anom, R_anom] = anom_meas(pos(1,k), pos(2,k));
    H_full           = [H; H_anom, zeros(2, 6)];
    R_full           = blkdiag(R, R_anom);

    Pdot = F * P + P * F' + Q - P * H_full' * (R_full \ (H_full * P));

    M              = L \ Pdot / L';
    S              = tril(M, -1) + 0.5 * diag(diag(M));
    Ldot           = L * S;
    Ldot_vec(:, k) = utils.L_to_lvec(Ldot);
end
end
