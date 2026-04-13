function [w] = compute_quadrature_weights(N)
% Compute Gauss-Lobatto quadrature weights.
%
% Inputs:
%   N - number of collocation intervals
%
% Outputs:
%   w - quadrature weights (1 x N+1)

[L, ~] = lpm.compute_legendre_polynomial(N);
tau    = lpm.compute_collocation_points(N);
P      = polyval(L, tau);
w      = 2 ./ (N * (N + 1) * (P .^ 2));
end