function lambda = lambda_from_filter(dt, T_filter)
%LAMBDA_FROM_FILTER  Decay factor equivalent to a first-order IIR low-pass.
%
%   lambda = lambda_from_filter(dt, T_filter)
%
%   The discrete filter y(k) = (1-a) y(k-1) + a x(k), a = dt/(T_filter+dt),
%   is an EWMA with lambda = 1 - a = T_filter/(T_filter + dt).

    lambda = T_filter ./ (T_filter + dt);
end
