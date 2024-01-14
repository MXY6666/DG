function tc = createTestCase

% viscosity coefficient 
tc.nu = 1.e-5;

% exact solution and its derivative
tc.u = @(x, t)0.5 * (1 - tanh(0.5 / sqrt(2 * tc.nu) * (x - 3 * sqrt(0.5 * tc.nu) * t)));

tc.u_x = @(x, t)-0.25 / sqrt(0.5 * tc.nu) * (1 - tanh(0.5 / sqrt(2 * tc.nu) * (x - 3 * sqrt(0.5 * tc.nu) * t)).^2);

% initial solution
tc.u0 = @(x)tc.u(x, 0);

% domain of computation
tc.dm = [-0.5, 2];


end