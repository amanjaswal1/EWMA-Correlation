function r = bond_return(y_prev, y_now, maturity, freq)
%BOND_RETURN  Monthly total return of a constant-maturity par bond.
%
%   r = bond_return(y_prev, y_now)
%   r = bond_return(y_prev, y_now, maturity, freq)
%
%   y_prev, y_now  Yields at the start and end of the month, as decimals
%                  (0.04 = 4 %). Vectors allowed.
%   maturity       Years (default 10).   freq  Coupons per year (default 2).
%
%   A par bond bought at y_prev pays coupon c = y_prev. One month later it
%   is repriced at y_now:
%       r = c/12 + P(c, y_now) - 1,
%       P(c, y) = (c/y)(1 - (1+y/f)^-n) + (1+y/f)^-n,   n = f * maturity.
%   Approximation: ignores the one month of roll-down in maturity.

    if nargin < 3, maturity = 10; end
    if nargin < 4, freq = 2; end
    c = y_prev;
    n = freq * maturity;
    i = y_now / freq;
    disc = (1 + i) .^ (-n);
    P = zeros(size(y_now));
    small = abs(i) < 1e-12;
    P(~small) = (c(~small) ./ y_now(~small)) .* (1 - disc(~small)) + disc(~small);
    P(small)  = c(small) * maturity + 1;        % limit y -> 0
    r = c / 12 + P - 1;
end
