function st = ewma_update(st, x)
%EWMA_UPDATE  Add one observation (1 x N) to a streaming EWMA estimator.
%
%   st = ewma_update(st, x)
%   st.mu is the current mean, st.S the current covariance matrix and
%   st.S(1,2)/sqrt(st.S(1,1)*st.S(2,2)) the current correlation.

    x = x(:).';
    st.k = st.k + 1;
    if st.warm
        a = max(1 - st.lambda, 1 / st.k);
    elseif st.k == 1
        a = 1;                                  % seed with the first value
    else
        a = 1 - st.lambda;
    end
    if st.demean
        st.mu = (1 - a) * st.mu + a * x;
    end
    d = x - st.mu;
    if ~st.warm && st.k == 1
        st.S = d.' * d;                         % = 0 when demeaning
    else
        st.S = (1 - a) * st.S + a * (d.' * d);
    end
end
