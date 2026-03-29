function Ldot_vec = cov_dynamics(psi, L_vec, Q)
% Compute time derivative of the Cholesky factor of the INS error covariance.
%
% The error state is [dp_N, dp_E, dv_N, dv_E, dpsi, bax, bay, bgz], where
% dp is position error [m], dv is velocity error [m/s], dpsi is heading
% error [rad], and bax/bay/bgz are accelerometer/gyroscope biases [m/s^2, rad/s].
%
% Inputs:
%   psi   - vehicle yaw (1 x n) [rad]
%   L_vec - lower Cholesky factor of error covariance, vectorized (36 x n)
%   Q     - process noise covariance matrix (8 x 8)
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

    Pdot = F * P + P * F' + Q;

    M              = L \ Pdot / L';
    S              = tril(M, -1) + 0.5 * diag(diag(M));
    Ldot           = L * S;
    Ldot_vec(:, k) = utils.L_to_lvec(Ldot);
end
end
