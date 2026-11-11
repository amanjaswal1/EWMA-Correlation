function p = ewma_params(lambda)
%EWMA_PARAMS  Interpretable numbers for an EWMA decay factor.
%
%   p = ewma_params(lambda)
%
%   p.halflife   periods until an observation's weight halves: log(0.5)/log(lambda)
%   p.mean_age   average age of the weighted data: lambda/(1-lambda)
%   p.n_eff      effective sample size (Kish): (1+lambda)/(1-lambda); the
%                sampling variance of an EWMA correlation is roughly
%                (1-rho^2)^2 / n_eff for normal data
%
%   For a first-order low-pass filter with time constant T sampled every dt,
%   lambda = T/(T + dt) (see LAMBDA_FROM_FILTER).

    p.lambda   = lambda;
    p.halflife = log(0.5) ./ log(lambda);
    p.mean_age = lambda ./ (1 - lambda);
    p.n_eff    = (1 + lambda) ./ (1 - lambda);
end
