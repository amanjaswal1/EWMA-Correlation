function test_warm_start()
% With warm start, the first K = 1/(1-lambda) estimates are exactly the
% cumulative sample means; afterwards the weight is the fixed 1 - lambda.
    rng(2);
    x = randn(100, 1);
    y = ewma_filter(x, 0.9, true);                % K = 10
    for t = 1:10
        assert(abs(y(t) - mean(x(1:t))) < 1e-14, sprintf('t=%d not the sample mean', t));
    end
    assert(abs(y(11) - (0.9 * y(10) + 0.1 * x(11))) < 1e-14, 'should switch to EWMA');
end
