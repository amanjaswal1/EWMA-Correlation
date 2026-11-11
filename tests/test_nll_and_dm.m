function test_nll_and_dm()
% GAUSS_NLL must match the textbook formula; DM must detect a clearly worse
% forecast and not reject when two forecasts are identical.
    rng(4);
    Sig = [1 0.5; 0.5 2];
    X = randn(2000, 2) * chol(Sig);
    H = repmat(Sig, [1 1 2000]);
    nll = gauss_nll(X, H);
    e = X(7, :).';
    ref = 0.5 * (2*log(2*pi) + log(det(Sig)) + e.' * (Sig \ e));
    assert(abs(nll(7) - ref) < 1e-12, 'NLL formula');
    nll_bad = gauss_nll(X, repmat(4 * Sig, [1 1 2000]));
    [t, p] = dm_test(nll, nll_bad);
    assert(t < -5 && p < 1e-6, 'DM should prefer the true covariance');
    [t0, p0] = dm_test(nll, nll);
    assert(isnan(t0) || (abs(t0) < 1e-12) || p0 == 1, 'identical losses');
end
