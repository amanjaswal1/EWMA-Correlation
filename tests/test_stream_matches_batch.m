function test_stream_matches_batch()
% Streaming EWMA_UPDATE must reproduce the vectorised EWMA_COV exactly,
% for every combination of demean / warm start.
    rng(1);
    X = randn(400, 3) * [1 0.3 0; 0 1 0.2; 0 0 1];
    for demean = [true false]
        for warm = [true false]
            S = ewma_cov(X, 0.9, demean, warm);
            st = ewma_init(3, 0.9, demean, warm);
            err = 0;
            for t = 1:size(X, 1)
                st = ewma_update(st, X(t, :));
                err = max(err, max(max(abs(st.S - S(:, :, t)))));
            end
            assert(err < 1e-12, sprintf('demean=%d warm=%d: max diff %.3g', demean, warm, err));
        end
    end
end
