function y = ewma_filter(x, lambda, warm)
%EWMA_FILTER  Exponentially weighted moving average of each column of x.
%
%   y = ewma_filter(x, lambda)
%   y = ewma_filter(x, lambda, warm)
%
%   y(t) = lambda * y(t-1) + (1 - lambda) * x(t)
%
%   warm = true (default): the weight is a_t = max(1 - lambda, 1/t), so the
%          first ~1/(1-lambda) values are plain cumulative means and the
%          estimate never depends on an arbitrary starting value.
%   warm = false: classic recursion started from y(1) = x(1).
%
%   This is the same first-order IIR low-pass filter used for the running
%   averages in a correlation-based phase-shift estimator; there
%   lambda = T_filter / (T_filter + dt).

    if nargin < 3, warm = true; end
    if lambda <= 0 || lambda >= 1
        error('ewma_filter:lambda', 'lambda must be in (0, 1).');
    end
    a = 1 - lambda;
    if isrow(x), x = x(:); flip_back = true; else, flip_back = false; end
    T = size(x, 1);
    y = zeros(size(x));
    if ~warm
        y = filter(a, [1, -lambda], x, lambda * x(1, :));
    else
        K = floor(1 / a);                      % last t with 1/t > a
        if K >= 1 && 1 / K <= a, K = K - 1; end
        if 1 / (K + 1) > a,      K = K + 1; end
        K = min(K, T);
        if K > 0
            y(1:K, :) = cumsum(x(1:K, :), 1) ./ (1:K)';
        end
        if K < T
            if K > 0, y_K = y(K, :); else, y_K = zeros(1, size(x, 2)); end
            y(K+1:T, :) = filter(a, [1, -lambda], x(K+1:T, :), lambda * y_K);
        end
    end
    if flip_back, y = y.'; end
end
