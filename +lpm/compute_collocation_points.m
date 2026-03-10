function [col] = compute_collocation_points(N)
% Compute Legendre-Gauss-Lobatto collocation points on [-1, 1].
%
% Inputs:
%   N - number of collocation intervals
%
% Outputs:
%   col - collocation points (1 x N+1)

% compute derivative of legendre polynomial
[~, Ldot] = lpm.compute_legendre_polynomial(N);

% form collocation points
col = reshape(sort([-1; roots(Ldot); 1]), 1, []);
end