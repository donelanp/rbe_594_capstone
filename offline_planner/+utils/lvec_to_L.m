function [L] = lvec_to_L(l_vec)
% Reconstruct a lower triangular matrix from a column vector.
%
% Inputs:
%   l_vec - lower triangular elements, column-major order (n*(n+1)/2 x 1)
%
% Outputs:
%   L - lower triangular matrix (n x n)

n                = (-1 + sqrt(1 + 8*numel(l_vec))) / 2;
L                = zeros(n);
L(tril(true(n))) = l_vec;
end
