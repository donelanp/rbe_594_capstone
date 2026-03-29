function L = lvec_to_L(l_vec)
% Reconstruct a lower triangular matrix from a column vector.
%
% Inputs:
%   l_vec - lower triangular elements, column-major order (36 x 1)
%
% Outputs:
%   L - lower triangular matrix (8 x 8)

L                  = zeros(8,8);
L(tril(true(8,8))) = l_vec;
end
