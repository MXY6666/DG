function tc = createTestCase(beta)

% advection coefficient
tc.beta = beta;

% initial solution
tc.u0 = @(x, y)(1 + sin(2 * pi * (x + y))) / 10;

% domain of computation
tc.dm = [0, 1, 0, 1];

% exact solution
tc.u = @(x, y, t)tc.u0(x - beta(1) * t, y - beta(2) * t);

end