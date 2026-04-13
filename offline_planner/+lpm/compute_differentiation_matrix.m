function [D] = compute_differentiation_matrix(N)
% Compute Legendre pseudospectral differentiation matrix.
%
% Inputs:
%   N - number of collocation intervals
%
% Outputs:
%   D - differentiation matrix (N+1 x N+1)

% compute legendre polynomial
[L, ~] = lpm.compute_legendre_polynomial(N);

% compute collocation points
col = lpm.compute_collocation_points(N);

% form differentiation matrix
D = zeros(N+1, N+1);
for jj=1:(N+1)
    for kk=1:(N+1)
        if jj ~= kk
            D(jj,kk) = polyval(L, col(jj)) / polyval(L, col(kk)) / (col(jj) - col(kk));
        elseif jj == 1
            D(jj,kk) = -N * (N+1) / 4;
        elseif jj == N+1
            D(jj,kk) = N * (N+1) / 4;
        end
    end
end
end