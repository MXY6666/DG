function tc = createTestCase(kappa, T, R, a, b)

% constant capillary coefficient
tc.kappa = kappa;

% temperature
tc.T = T;

% universal gas constant
tc.R = R;

% positive constants depending on the type of fluid
tc.a = a;
tc.b = b;

% initial solutions
tc.rho0 = @(x)0.5 + 0.2 * sin(x);
tc.u0   = @(x)cos(x);

% exact solutions and their derivatives
tc.rho = @(x, t)0.5 + 0.2 * exp(-t) * sin(x);
tc.u   = @(x, t)exp(-t) * cos(x);

tc.rho_t   = @(x, t)-0.2 * exp(-t) * sin(x);
tc.rho_x   = @(x, t) 0.2 * exp(-t) * cos(x);
tc.rho_xxx = @(x, t)-0.2 * exp(-t) * cos(x);

tc.u_t = @(x, t)-exp(-t) * cos(x);
tc.u_x = @(x, t)-exp(-t) * sin(x);

% domain of computation
tc.dm = [0, 2 * pi];

end