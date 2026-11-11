function st = ewma_init(N, lambda, demean, warm)
%EWMA_INIT  State for a streaming (one observation at a time) EWMA estimator.
%
%   st = ewma_init(N, lambda)
%   st = ewma_init(N, lambda, demean, warm)
%
%   Use with EWMA_UPDATE in a live loop; gives the same numbers as EWMA_COV.

    if nargin < 3 || isempty(demean), demean = true; end
    if nargin < 4 || isempty(warm),   warm = true;   end
    st = struct('N', N, 'lambda', lambda, 'demean', demean, 'warm', warm, ...
                'k', 0, 'mu', zeros(1, N), 'S', zeros(N, N));
end
