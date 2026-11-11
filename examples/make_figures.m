% MAKE_FIGURES  Regenerate every figure and number quoted in README.md.
%
%   run('examples/make_figures.m')  will save -->  results/figures/*.png + printed numbers

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
out_dir = fullfile(here, '..', 'results', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Real data ---------------------------------------------------------------
d = load_shiller();
X = d.X; yr = d.year; T = size(X, 1);
fprintf('\nData: %d monthly returns, %s to %s\n', T, d.date{1}, d.date{end});
fprintf('  annualised vol: stocks %.1f %%, 10y bonds %.1f %%\n', 100 * std(X) * sqrt(12));

% Walk-forward selection of lambda (scored from 1963, tested from 1990)
lam_grid = [0.80:0.01:0.99, 0.995];
eval_start = 1963; test_start = 1990;
wf = walk_forward_lambda(X, yr, lam_grid, eval_start, test_start);
train = wf.eval & ~wf.test;
[~, j_train] = min(mean(wf.nll(train, :), 1));
lam_star = lam_grid(j_train);
p_star = ewma_params(lam_star);
fprintf('\nWalk-forward lambda selection (Gaussian NLL of one-step forecasts):\n');
fprintf('  best on %d-%d: lambda = %.2f (half-life %.1f months, n_eff %.0f)\n', ...
        eval_start, test_start - 1, lam_star, p_star.halflife, p_star.n_eff);
fprintf('  yearly walk-forward choices %d-%d: %.2f .. %.2f\n', ...
        wf.test_years(1), wf.test_years(end), min(wf.lam_path), max(wf.lam_path));
[~, j_hind] = min(mean(wf.nll(wf.test, :), 1));
fprintf('  (hindsight only) best fixed lambda on the test period: %.2f\n', lam_grid(j_hind));

% Benchmarks on the same test months
models = {'EWMA, walk-forward lambda', sprintf('EWMA, fixed lambda = %.2f (chosen pre-1990)', lam_star), ...
          'EWMA, RiskMetrics lambda = 0.97', 'rolling 60-month window', ...
          'rolling 120-month window', 'expanding window (all past data)'};
loss = NaN(T, numel(models));
loss(:, 1) = wf.nll_adapt;
loss(:, 2) = wf.nll(:, j_train);
loss(:, 3) = wf.nll(:, abs(lam_grid - 0.97) < 1e-12);
wins = [60, 120, Inf];
for k = 1:3
    [S, mu] = rolling_cov(X, wins(k));
    [H, M] = one_step_ahead(S, mu);
    loss(:, 3 + k) = gauss_nll(X, H, M);
end
te = wf.test;
fprintf('\nOut-of-sample %d-%d (%d months), mean NLL per month (lower = better):\n', ...
        test_start, floor(yr(end)), sum(te));
fprintf('  %-48s %9s %9s %9s\n', 'model', 'mean NLL', 'DM t', 'p-value');
for k = 1:numel(models)
    if k == 1
        fprintf('  %-48s %9.4f %9s %9s\n', models{k}, mean(loss(te, k)), '-', '-');
    else
        [t_dm, p_dm] = dm_test(loss(te, 1), loss(te, k));
        fprintf('  %-48s %9.4f %9.2f %9.3f\n', models{k}, mean(loss(te, k)), t_dm, p_dm);
    end
end
fprintf('  (DM t < 0: walk-forward EWMA has lower loss than that model)\n');

%% Figure 1 - stock-bond correlation since 1958 -----------------------------
S = ewma_cov(X, lam_star);
rho_ewma = ewma_corr(S);
S60 = rolling_cov(X, 60);
rho_60 = ewma_corr(S60);
show = yr >= 1958;

eras = [1958 1970; 1970 2000; 2000 2022; 2022 yr(end) + 1/12];
fprintf('\nSample correlation of monthly stock and bond returns by era:\n');
era_rho = zeros(size(eras, 1), 1);
for e = 1:size(eras, 1)
    idx = yr >= eras(e, 1) & yr < eras(e, 2);
    c = corrcoef(X(idx, 1), X(idx, 2));
    era_rho(e) = c(1, 2);
    fprintf('  %d-%d: %+.2f  (%d months)\n', eras(e, 1), floor(min(eras(e, 2) - 1/12, yr(end))), era_rho(e), sum(idx));
end

fig = figure('Visible', 'off', 'Color', 'w'); hold on;
shade = [0.93 0.95 0.98; 0.98 0.95 0.92];
for e = 1:size(eras, 1)
    fill([eras(e, 1) eras(e, 2) eras(e, 2) eras(e, 1)], [-1 -1 1 1], shade(1 + mod(e, 2), :), ...
         'EdgeColor', 'none');
    plot(eras(e, :), era_rho(e) * [1 1], '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5);
end
plot([1958 2024], [0 0], 'k-', 'LineWidth', 0.5);
h1 = plot(yr(show), rho_60(show), '-', 'Color', palette('orange'), 'LineWidth', 1.3);
h2 = plot(yr(show), rho_ewma(show), '-', 'Color', palette('blue'), 'LineWidth', 1.8);
h3 = plot(NaN, NaN, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5);
grid on; box on; xlim([1958 2024]); ylim([-0.8 0.8]);
set(gca, 'Layer', 'top');
xlabel('year'); ylabel('correlation');
legend([h2 h1 h3], {sprintf('EWMA, \\lambda = %.2f (half-life %.0f months)', lam_star, p_star.halflife), ...
        'rolling 60-month window', 'era average (sample correlation)'}, ...
       'Location', 'southoutside', 'NumColumns', 3);
title('U.S. stocks vs 10-year Treasuries: monthly return correlation regimes since 1958');
save_png(fig, fullfile(out_dir, 'fig1_stock_bond_correlation.png'), 1100, 600);

%% Figure 2 - known-truth simulation: the bias-variance trade-off ----------
rng(1);
A = 0.6; P = 500; Tsim = 6000; burn = 1000; M = 100;
rho_true = A * sin(2*pi*(1:Tsim)'/P);
lams = 1 - logspace(log10(0.25), log10(0.003), 25);
mse_mc = zeros(M, numel(lams));
for m = 1:M
    Xs = sim_corr_returns(rho_true);
    for j = 1:numel(lams)
        r = ewma_corr(ewma_cov(Xs, lams(j)));
        mse_mc(m, j) = mean((r(burn+1:end) - rho_true(burn+1:end)) .^ 2);
    end
end
[mse_th, b2_th, var_th] = ewma_corr_mse_theory(lams, A, P);
[~, jm] = min(mean(mse_mc, 1));
lam_fine = 1 - logspace(log10(0.25), log10(0.003), 2000);
[~, jt] = min(ewma_corr_mse_theory(lam_fine, A, P));
fprintf('\nSimulation (true rho = %.1f sin(2 pi t / %d), %d paths x %d obs):\n', A, P, M, Tsim);
fprintf('  MSE-optimal lambda: Monte Carlo %.3f, theory %.3f\n', lams(jm), lam_fine(jt));
fprintf('  max |MC/theory - 1| over the lambda grid: %.1f %%\n', 100 * max(abs(mean(mse_mc, 1) ./ mse_th - 1)));

Xs = sim_corr_returns(rho_true);
fig = figure('Visible', 'off', 'Color', 'w');
subplot(2, 1, 1); hold on;
win = 3001:4000;
plot(win, rho_true(win), 'k-', 'LineWidth', 2.5);
show_lams = [0.90, lams(jm), 0.99];
cols = {palette('orange'), palette('blue'), palette('purple')};
for j = 1:3
    r = ewma_corr(ewma_cov(Xs, show_lams(j)));
    plot(win, r(win), '-', 'Color', cols{j}, 'LineWidth', 1.3);
end
grid on; box on; ylim([-1 1]); xlim([win(1) win(end)]);
xlabel('time [observations]'); ylabel('correlation');
legend({'true correlation', '\lambda = 0.90: noisy', sprintf('\\lambda = %.3f: MSE-optimal', lams(jm)), ...
        '\lambda = 0.99: smooth but late and damped'}, 'Location', 'southoutside', 'NumColumns', 2);
title('(a) One simulated path');

subplot(2, 1, 2);
x = 1 - lams;
loglog(x, mean(mse_mc, 1), 'o', 'Color', palette('blue'), 'MarkerFaceColor', palette('blue'), 'MarkerSize', 5); hold on;
loglog(x, mse_th, '-', 'Color', 'k', 'LineWidth', 1.6);
loglog(x, b2_th, '--', 'Color', palette('red'), 'LineWidth', 1.3);
loglog(x, var_th, '--', 'Color', palette('green'), 'LineWidth', 1.3);
grid on; box on; xlim([min(x) max(x)]); ylim([1e-3 0.3]);
set(gca, 'XDir', 'reverse');
xlabel('1 - \lambda   (right = shorter memory)'); ylabel('mean squared error');
legend({'Monte Carlo (100 paths)', 'theory: bias^2 + variance', 'theory: bias^2 (lag)', ...
        'theory: variance (noise)'}, 'Location', 'southoutside', 'NumColumns', 2);
title('(b) Bias-variance trade-off: closed-form theory vs simulation');
save_png(fig, fullfile(out_dir, 'fig2_bias_variance.png'), 1000, 1000);

%% Figure 3 - lag after a correlation break --------------------------------
rng(2);
Tb = 600; t0 = 300; rho_old = 0.5; rho_new = -0.3; Mb = 400;
rho_step = rho_old * ones(Tb, 1); rho_step(t0+1:end) = rho_new;
step_lams = [0.90, 0.97];
fig = figure('Visible', 'off', 'Color', 'w'); hold on;
cols = {palette('orange'), palette('blue')};
hl = [];
for j = 1:2
    R = zeros(Tb, Mb);
    for m = 1:Mb
        R(:, m) = ewma_corr(ewma_cov(sim_corr_returns(rho_step), step_lams(j)));
    end
    mR = mean(R, 2); sR = std(R, 0, 2);
    k = (t0 - 50):Tb;
    fill([k fliplr(k)], [mR(k) + sR(k); flipud(mR(k) - sR(k))]', cols{j}, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');
    hl(end+1) = plot(k, mR(k), '-', 'Color', cols{j}, 'LineWidth', 2); %#ok<SAGROW>
    kk = t0:Tb;
    hl(end+1) = plot(kk, rho_new + (rho_old - rho_new) * step_lams(j) .^ (kk - t0), ':', ...
                     'Color', 'k', 'LineWidth', 1.5); %#ok<SAGROW>
    p = ewma_params(step_lams(j));
    fprintf('  step response lambda = %.2f: half-life %.1f obs; MC mean crosses the midpoint after %d obs\n', ...
            step_lams(j), p.halflife, find(mR(t0+1:end) < (rho_old + rho_new)/2, 1));
end
plot([t0 t0], [-1 1], 'k--');
grid on; box on; ylim([-0.8 0.9]); xlim([t0 - 50, Tb]);
xlabel('time [observations]'); ylabel('estimated correlation');
legend(hl([1 3 2]), {'\lambda = 0.90: mean of 400 paths (band = \pm1 sd)', ...
       '\lambda = 0.97: mean of 400 paths (band = \pm1 sd)', ...
       'theory: \rho_{new} + (\rho_{old} - \rho_{new}) \lambda^{t - t_0}'}, 'Location', 'northeast');
title('Correlation break 0.5 \rightarrow -0.3: the estimate decays like \lambda^t');
save_png(fig, fullfile(out_dir, 'fig3_step_response.png'), 1000, 560);

%% Figure 4 - choosing lambda on real data, out of sample -------------------
fig = figure('Visible', 'off', 'Color', 'w');
subplot(1, 2, 1); hold on;
nll_tr = mean(wf.nll(train, :), 1);
nll_te = mean(wf.nll(te, :), 1);
y_lo = min([nll_tr, nll_te, mean(loss(te, 4:6), 1)]) - 0.02;
y_hi = y_lo + 0.30;
plot(lam_grid, nll_tr, 'o-', 'Color', palette('grey'), 'LineWidth', 1.5, 'MarkerSize', 4);
plot(lam_grid, nll_te, 'o-', 'Color', palette('blue'), 'LineWidth', 1.5, 'MarkerSize', 4);
plot(lam_star * [1 1], [y_lo y_hi], '--', 'Color', palette('grey'));
plot(0.97 * [1 1], [y_lo y_hi], ':', 'Color', 'k', 'LineWidth', 1.2);
bench = {'rolling 60 m (test)', 'rolling 120 m (test)', 'expanding (test)'};
for k = 4:6
    v = mean(loss(te, k));
    plot([0.80 0.87], v * [1 1], '-', 'Color', palette('orange'), 'LineWidth', 1.2);
    text(0.802, v + 0.009, bench{k - 3}, 'Color', palette('orange'), 'FontSize', 9);
end
j86 = find(abs(lam_grid - 0.86) < 1e-12);
j90 = find(abs(lam_grid - 0.90) < 1e-12);
text(0.856, nll_tr(j86), 'fit 1963-89', 'Color', palette('grey'), 'FontSize', 10, ...
     'HorizontalAlignment', 'right', 'FontWeight', 'bold');
text(0.906, nll_te(j90), 'test 1990-2023', 'Color', palette('blue'), 'FontSize', 10, ...
     'FontWeight', 'bold');
text(lam_star - 0.004, y_hi - 0.015, sprintf('chosen %.2f', lam_star), 'HorizontalAlignment', 'right', 'FontSize', 9);
text(0.973, y_hi - 0.015, 'RiskMetrics 0.97', 'FontSize', 9, 'HorizontalAlignment', 'right', ...
     'Rotation', 90, 'VerticalAlignment', 'top');
grid on; box on; xlim([0.80 1.0]); ylim([y_lo y_hi]);
xlabel('\lambda'); ylabel('mean NLL per month (lower = better)');
title('(a) Loss vs \lambda');
subplot(1, 2, 2);
stairs([wf.test_years; wf.test_years(end) + 1], [wf.lam_path; wf.lam_path(end)], ...
       '-', 'Color', palette('blue'), 'LineWidth', 2);
grid on; box on; ylim([0.90 1.0]); xlim([wf.test_years(1), wf.test_years(end) + 1]);
xlabel('year'); ylabel('\lambda used that year');
title('(b) Walk-forward choice, refit each January');
save_png(fig, fullfile(out_dir, 'fig4_walk_forward.png'), 1100, 480);
fprintf('\n');
