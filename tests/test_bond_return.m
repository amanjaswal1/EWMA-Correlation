function test_bond_return()
% Par bond: unchanged yield earns exactly the coupon; a 1 bp rise loses
% about modified duration x 1 bp; zero-yield limit is handled.
    assert(abs(bond_return(0.04, 0.04) - 0.04/12) < 1e-15, 'unchanged yield');
    D_mod = (1 - 1.02^-20) / 0.04;            % par bond: (1/y)(1-(1+y/2)^-2n) = 8.18 y
    r = bond_return(0.04, 0.0401) - 0.04/12;
    assert(abs(r / -1e-4 - D_mod) < 0.01, sprintf('duration %.3f vs %.3f', r / -1e-4, D_mod));
    assert(isfinite(bond_return(0.01, 0)), 'zero yield');
end
