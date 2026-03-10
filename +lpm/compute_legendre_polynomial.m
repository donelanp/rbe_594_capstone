function [L, Ldot] = compute_legendre_polynomial(N)
% Compute Nth-order Legendre polynomial and its derivative.
%
% Inputs:
%   N - polynomial order
%
% Outputs:
%   L    - Legendre polynomial coefficients (1 x N+1)
%   Ldot - derivative polynomial coefficients (1 x N)

% populate output as array of coefficients
L = zeros(1, N + 1);

for k=0:floor(N / 2)
    power          = N - 2*k;
    coeff          = 2^(-N) * (-1)^k * nchoosek(N, k) * nchoosek(2*N - 2*k, N);
    L(end - power) = L(end - power) + coeff;
end

Ldot = polyder(L);
end