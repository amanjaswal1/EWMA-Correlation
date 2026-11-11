function [S, mu] = ewma_cov(X, lambda, demean, warm)
%EWMA_COV  Exponentially weighted covariance matrices through time.
%
%   [S, mu] = ewma_cov(X, lambda)
%   [S, mu] = ewma_cov(X, lambda, demean, warm)
%
%   X        T x N returns (one column per asset).
%   lambda   Decay factor in (0, 1); RiskMetrics uses 0.94 (daily) and
%            0.97 (monthly).
%   demean   true (default): deviations from the EWMA mean, as in the
%            phase-shift estimator. false: zero-mean (RiskMetrics) form.
%   warm     true (default): cumulative-mean start-up, see EWMA_FILTER.
%
%   S        N x N x T. S(:,:,t) uses data up to and including t, so the
%            one-step-ahead forecast for t+1 is S(:,:,t).
%   mu       T x N EWMA means (zeros if demean = false).
%
%   Every entry of S is an EWMA of a cross-product series, so the whole
%   history is computed with N(N+1)/2 vectorised filter calls.

    if nargin < 3 || isempty(demean), demean = true; end
    if nargin < 4 || isempty(warm),   warm = true;   end
    [T, N] = size(X);
    if demean
        mu = ewma_filter(X, lambda, warm);
    else
        mu = zeros(T, N);
    end
    D = X - mu;
    S = zeros(N, N, T);
    for i = 1:N
        for j = i:N
            s = ewma_filter(D(:, i) .* D(:, j), lambda, warm);
            S(i, j, :) = reshape(s, 1, 1, T);
            S(j, i, :) = reshape(s, 1, 1, T);
        end
    end
end
