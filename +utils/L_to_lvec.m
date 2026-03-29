function l_vec = L_to_lvec(L)
% Extract lower-triangular elements of a matrix into a column vector.
%
% Inputs:
%   L - lower triangular matrix (8 x 8)
%
% Outputs:
%   l_vec - lower triangular elements, column-major order (36 x 1)

l_vec = L(tril(true(8,8)));
end
