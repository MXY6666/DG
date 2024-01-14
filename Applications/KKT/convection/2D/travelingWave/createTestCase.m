function tc = createTestCase(beta, KKT)

% advection coefficient
tc.beta = beta;

% initial solution
tc.u0 = @(x, y)max(cos(2 * pi / 3 * x) .* cos(2 * pi / 3 * y), KKT.Umin);

% domain of computation
tc.dm = [0, 3, 0, 3];

% exact solution
tc.u = @(x, y, t)tc.u0(x - beta(1) * t, y - beta(2) * t);

end