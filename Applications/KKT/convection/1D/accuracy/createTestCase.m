function tc = createTestCase(beta)

% advection coefficient
tc.beta = beta;

% initial solution
tc.u0 = @(x)(1 + sin(2 * pi * x)) / 10;

% domain of computation
tc.dm = [0, 1];

% exact solution
tc.u = @(x, t)tc.u0(x - beta * t);

end