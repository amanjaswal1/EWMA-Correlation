function test_mse_theory()
% The bias-variance formula must match Monte Carlo within 10 % for a
% sinusoidal true correlation (A = 0.6, period 500).
    rng(3);
    A = 0.6; P = 500; T = 5000; burn = 1000; M = 20;
    rho = A * sin(2*pi*(1:T)'/P);
    for lam = [0.9 0.95 0.98]
        mse = zeros(M, 1);
        for m = 1:M
            r = ewma_corr(ewma_cov(sim_corr_returns(rho), lam));
            mse(m) = mean((r(burn+1:end) - rho(burn+1:end)) .^ 2);
        end
        th = ewma_corr_mse_theory(lam, A, P);
        assert(abs(mean(mse) / th - 1) < 0.10, sprintf('lambda=%.2f: MC %.4f vs theory %.4f', lam, mean(mse), th));
    end
end
