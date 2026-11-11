function nll = gauss_nll(X, H, mu)
%GAUSS_NLL  Per-period Gaussian negative log-likelihood of covariance forecasts.
%
%   nll = gauss_nll(X, H)
%   nll = gauss_nll(X, H, mu)
%
%   X    T x N realised returns.
%   H    N x N x T forecast covariance for each period (made BEFORE seeing X(t,:)).
%   mu   T x N forecast means (default zeros).
%   nll  T x 1: 0.5*(N*log(2*pi) + log det H_t + e_t' inv(H_t) e_t).
%
%   Lower is genrally better. It is a proper scoring rule for covariance forecasts:
%   in expectation it is minimised by the true covariance matrix.

    [T, N] = size(X);
    if nargin < 3 || isempty(mu), mu = zeros(T, N); end
    nll = NaN(T, 1);
    for t = 1:T
        Ht = H(:, :, t);
        if any(~isfinite(Ht(:))), continue, end
        e = (X(t, :) - mu(t, :)).';
        [R, p] = chol(Ht);
        if p > 0, continue, end                    % not positive definite
        z = R.' \ e;
        nll(t) = 0.5 * (N * log(2*pi) + 2 * sum(log(diag(R))) + z.' * z);
    end
end
