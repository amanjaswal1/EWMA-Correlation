function test_params()
% Known values: RiskMetrics daily lambda = 0.94 has a half-life of 11.2
% days; the phase estimator's filter (dt = 1e-4 s, T = 0.1 s) is an EWMA
% with lambda = 0.999 and a half-life of 69.3 ms.
    p = ewma_params(0.94);
    assert(abs(p.halflife - 11.2023) < 1e-3, 'half-life of 0.94');
    assert(abs(p.n_eff - 1.94 / 0.06) < 1e-12, 'n_eff');
    lam = lambda_from_filter(1e-4, 0.1);
    assert(abs(lam - 0.1 / 0.1001) < 1e-15, 'filter mapping');
    assert(abs(ewma_params(lam).halflife * 1e-4 - 0.0693) < 1e-4, 'filter half-life in s');
end
