function X = sim_corr_returns(rho, sigma)
%SIM_CORR_RETURNS  Bivariate normal returns with a given correlation path.
%
%   X = sim_corr_returns(rho)
%   X = sim_corr_returns(rho, sigma)
%
%   rho    T x 1 true correlation at each time (|rho| < 1).
%   sigma  1 x 2 volatilities (default [1 1]).
%   X      T x 2 returns, zero mean: x1 = s1*z1,
%          x2 = s2*(rho*z1 + sqrt(1-rho^2)*z2), z ~ iid N(0,1).
%   Uses the global random stream; call rng(seed) first for reproducibility.

    if nargin < 2, sigma = [1 1]; end
    rho = rho(:);
    Z = randn(numel(rho), 2);
    X = [sigma(1) * Z(:, 1), ...
         sigma(2) * (rho .* Z(:, 1) + sqrt(1 - rho .^ 2) .* Z(:, 2))];
end
