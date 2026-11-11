function wf = walk_forward_lambda(X, year, lam_grid, eval_start, test_start)
%WALK_FORWARD_LAMBDA  Choose the EWMA decay factor without look-ahead.
%
%   wf = walk_forward_lambda(X, year, lam_grid, eval_start, test_start)
%
%   X           T x N returns;  year  T x 1 decimal year of each return.
%   lam_grid    Candidate decay factors.
%   eval_start  First year scored (earlier data only warms up the filters).
%   test_start  First out-of-sample year.
%
%   Every candidate produces one-step-ahead covariance forecasts, scored by
%   the Gaussian negative log-likelihood (GAUSS_NLL). Then, at the start of
%   each test year, the decay factor with the lowest average loss over all
%   scored data BEFORE that year is selected and used for the next 12 months.
%
%   wf.nll         T x L loss of every fixed lambda
%   wf.nll_adapt   T x 1 loss of the walk-forward choice (NaN before test)
%   wf.test_years, wf.lam_path   lambda chosen for each test year
%   wf.eval, wf.test             logical masks

    L = numel(lam_grid);
    T = size(X, 1);
    wf.nll = NaN(T, L);
    for j = 1:L
        [S, mu] = ewma_cov(X, lam_grid(j));
        [H, M] = one_step_ahead(S, mu);
        wf.nll(:, j) = gauss_nll(X, H, M);
    end
    wf.eval = year >= eval_start;
    wf.test = year >= test_start;
    wf.test_years = unique(floor(year(wf.test)));
    wf.lam_path = zeros(size(wf.test_years));
    wf.nll_adapt = NaN(T, 1);
    for k = 1:numel(wf.test_years)
        y = wf.test_years(k);
        fit = wf.eval & year < y;
        [~, j] = min(mean(wf.nll(fit, :), 1));
        wf.lam_path(k) = lam_grid(j);
        idx = floor(year) == y;
        wf.nll_adapt(idx) = wf.nll(idx, j);
    end
    wf.lam_grid = lam_grid;
end
