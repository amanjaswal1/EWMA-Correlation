function [mse, bias2, variance] = ewma_corr_mse_theory(lambda, A, period)
%EWMA_CORR_MSE_THEORY  Approximate MSE of an EWMA correlation estimate when
%the true correlation oscillates as rho_t = A*sin(2*pi*t/period).
%
%   [mse, bias2, variance] = ewma_corr_mse_theory(lambda, A, period)
%
%   bias2     Tracking error of the filter. The EWMA of the cross-product is
%             a linear filter H(w) = (1-lambda)/(1 - lambda*exp(-1i*w)), so
%             the mean estimate is the truth passed through H and
%             mean(bias^2) = A^2/2 * |1 - H(w)|^2,  w = 2*pi/period.
%   variance  Sampling noise. For normal data the variance of a correlation
%             estimate is (1-rho^2)^2 / n_eff with n_eff = (1+lambda)/(1-lambda);
%             averaged over the cycle, E[(1-rho^2)^2] = 1 - A^2 + 3A^4/8.
%   mse       bias2 + variance.
%
%   Small lambda: little lag, lots of noise. Large lambda: smooth but late.

    w = 2 * pi / period;
    H = (1 - lambda) ./ (1 - lambda .* exp(-1i * w));
    bias2 = (A^2 / 2) * abs(1 - H) .^ 2;
    n_eff = (1 + lambda) ./ (1 - lambda);
    variance = (1 - A^2 + 3 * A^4 / 8) ./ n_eff;
    mse = bias2 + variance;
end
