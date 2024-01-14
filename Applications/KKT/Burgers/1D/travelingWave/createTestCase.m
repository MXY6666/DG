function tc = createTestCase(KKT)

% initial solution
tc.u0 = @(x)max(cos(pi * x), KKT.Umin);

% exact solution
tc.u = @(x, t)fzero(@(v)v - tc.u0(x - v .* t), 0);

% domain of computation
tc.dm = [-1, 1];


end