function [l_vec] = L_to_lvec(L)
% Extract lower-triangular elements of a matrix into a column vector.
%
% Inputs:
%   L - lower triangular matrix (n x n)
%
% Outputs:
%   l_vec - lower triangular elements, column-major order (n*(n+1)/2 x 1)

l_vec = L(tril(true(size(L))));
end
