function tc = createTestCase(beta, KKT)

% advection coefficient
tc.beta = beta;

% initial solution
tc.u0 = @(x)max(cos(2 * pi * x / 10), KKT.Umin);

% domain of computation
tc.dm = [0, 10];

% exact solution
tc.u = @(x, t)tc.u0(x - beta * t);

end