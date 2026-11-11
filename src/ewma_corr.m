function rho = ewma_corr(S, i, j)
%EWMA_CORR  Correlation series from a covariance history.
%
%   rho = ewma_corr(S)         % pair (1,2)
%   rho = ewma_corr(S, i, j)   % pair (i,j)
%
%   S is N x N x T (from EWMA_COV); rho is T x 1.

    if nargin < 2, i = 1; j = 2; end
    sij = squeeze(S(i, j, :));
    sii = squeeze(S(i, i, :));
    sjj = squeeze(S(j, j, :));
    rho = sij ./ sqrt(sii .* sjj);
    rho(sii <= 0 | sjj <= 0) = NaN;
end
