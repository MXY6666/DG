function tc = createTestCase

% initial solution
tc.u0 = @(x)1 + sin(2 * pi * x);

% exact solution
tc.u = @(x, t)fzero(@(v)v - tc.u0(x - v .* t), 0);

% domain of computation
tc.dm = [0, 1];


end