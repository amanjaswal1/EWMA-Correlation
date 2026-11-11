function [t_stat, p_value, mean_diff] = dm_test(loss_a, loss_b, lags)
%DM_TEST  Diebold-Mariano test of equal predictive accuracy.
%
%   [t_stat, p_value, mean_diff] = dm_test(loss_a, loss_b)
%   [t_stat, p_value, mean_diff] = dm_test(loss_a, loss_b, lags)
%
%   d_t = loss_a(t) - loss_b(t). H0: E[d_t] = 0. The variance of the mean
%   uses a Newey-West (Bartlett) estimator with LAGS lags (default
%   floor(4*(T/100)^(2/9))). Negative t_stat = model A has lower loss.
%   p_value is two-sided, from the normal approximation.

    d = loss_a(:) - loss_b(:);
    d = d(isfinite(d));
    T = numel(d);
    if nargin < 3 || isempty(lags)
        lags = floor(4 * (T / 100) ^ (2 / 9));
    end
    mean_diff = mean(d);
    u = d - mean_diff;
    lrv = sum(u .^ 2) / T;
    for L = 1:lags
        w = 1 - L / (lags + 1);
        lrv = lrv + 2 * w * sum(u(L+1:end) .* u(1:end-L)) / T;
    end
    t_stat = mean_diff / sqrt(lrv / T);
    p_value = erfc(abs(t_stat) / sqrt(2));          % 2 * (1 - Phi(|t|))
end
