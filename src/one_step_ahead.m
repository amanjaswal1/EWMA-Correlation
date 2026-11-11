function [H, M] = one_step_ahead(S, mu)
%ONE_STEP_AHEAD  Turn filtered estimates into forecasts for the next period.
%
%   [H, M] = one_step_ahead(S, mu)
%
%   S(:,:,t), mu(t,:) use data up to t. The forecast for period t is the
%   estimate available at t-1, so H(:,:,t) = S(:,:,t-1), M(t,:) = mu(t-1,:);
%   period 1 has no forecast (NaN). Keeps the evaluation free of look-ahead.

    [N, ~, T] = size(S);
    H = cat(3, NaN(N, N), S(:, :, 1:T-1));
    M = [NaN(1, N); mu(1:T-1, :)];
end
