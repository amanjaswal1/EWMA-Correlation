% DEMO_STREAMING  Live (one month at a time) EWMA correlation on the data.
%
%   Shows the streaming API that a trading or risk system would call each
%   period, and prints the latest stock-bond correlation and volatilities.

clear;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));

d = load_shiller();
lambda = 0.95;
st = ewma_init(2, lambda);
for t = 1:size(d.X, 1)
    st = ewma_update(st, d.X(t, :));        % one new month of returns
end

S = st.S;
rho = S(1, 2) / sqrt(S(1, 1) * S(2, 2));
p = ewma_params(lambda);
fprintf('lambda = %.2f (half-life %.1f months, effective sample %.0f months)\n', ...
        lambda, p.halflife, p.n_eff);
fprintf('as of %s: stock vol %.1f %%/yr, bond vol %.1f %%/yr, correlation %+.2f\n', ...
        d.date{end}, 100 * sqrt(12 * S(1, 1)), 100 * sqrt(12 * S(2, 2)), rho);
