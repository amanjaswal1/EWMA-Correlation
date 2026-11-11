function test_data()
% The shipped data loads, has the expected span, and gives sensible returns.
    d = load_shiller();
    assert(size(d.X, 1) == 842, sprintf('expected 842 returns, got %d', size(d.X, 1)));
    assert(strcmp(d.date{1}, '1953-05') && strcmp(d.date{end}, '2023-06'), 'date range');
    assert(all(isfinite(d.X(:))), 'non-finite returns');
    assert(all(diff(d.year) > 0), 'dates not increasing');
    vol = std(d.X) * sqrt(12);
    assert(vol(1) > 0.10 && vol(1) < 0.20 && vol(2) > 0.03 && vol(2) < 0.12, 'implausible volatility');
end
