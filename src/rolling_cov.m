function [S, mu] = rolling_cov(X, window)
%ROLLING_COV  Equal-weight sample covariance over the last WINDOW periods.
%
%   [S, mu] = rolling_cov(X, window)   % N x N x T and T x N, NaN until t >= window
%   [S, mu] = rolling_cov(X, Inf)      % expanding window (all data so far, t >= 2)
%
%   Benchmark for the EWMA estimator; S(:,:,t) and mu(t,:) use data up to
%   and including t.

    [T, N] = size(X);
    S = NaN(N, N, T);
    mu = NaN(T, N);
    for t = 1:T
        if isinf(window)
            if t < 2, continue, end
            Z = X(1:t, :);
        else
            if t < window, continue, end
            Z = X(t - window + 1:t, :);
        end
        S(:, :, t) = cov(Z);
        mu(t, :) = mean(Z, 1);
    end
end
