function test_sinusoid_phase()
% Bridge to the phase-shift estimator: for two sinusoids shifted by phi,
% the EWMA correlation equals cos(phi).
    t = (0:1e-4:1).';
    phi = 60 * pi / 180;
    X = [sin(2*pi*50*t), 0.3 + 2*sin(2*pi*50*t - phi)];
    S = ewma_cov(X, lambda_from_filter(1e-4, 0.1));
    rho = ewma_corr(S);
    est = mean(rho(end - 199:end));               % average over one cycle
    assert(abs(est - cos(phi)) < 2e-3, sprintf('rho = %.4f, cos(phi) = %.4f', est, cos(phi)));
end
