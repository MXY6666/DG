function tc = createTestCase(Re, We, T, R, a, b)

% Reynolds number
tc.Re = Re;

% Weber number
tc.We = We;

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
tc.m0   = @(x)tc.rho0(x) .* tc.u0(x);

% exact solutions and their derivatives
tc.rho = @(x, t)0.5 + 0.2 * exp(-t) * sin(x);
tc.u   = @(x, t)exp(-t) * cos(x);

tc.rho_t   = @(x, t)-0.2 * exp(-t) * sin(x);
tc.rho_x   = @(x, t) 0.2 * exp(-t) * cos(x);
tc.rho_xxx = @(x, t)-0.2 * exp(-t) * cos(x);

tc.u_t  = @(x, t)-exp(-t) * cos(x);
tc.u_x  = @(x, t)-exp(-t) * sin(x);
tc.u_xx = @(x, t)-exp(-t) * cos(x);

% % initial solutions
% tc.rho0 = @(x)(1.2 + sin(x)) / 10;
% tc.u0   = @(x)cos(x);
% tc.m0   = @(x)tc.rho0(x) .* tc.u0(x);
% 
% % exact solutions and their derivatives
% tc.rho = @(x, t)(1.2 + sin(x + t)) / 10;
% tc.u   = @(x, t)cos(x + t);
% 
% tc.rho_t   = @(x, t) cos(x + t) / 10;
% tc.rho_x   = @(x, t) cos(x + t) / 10;
% tc.rho_xxx = @(x, t)-cos(x + t) / 10;
% 
% tc.u_t  = @(x, t)-sin(x + t);
% tc.u_x  = @(x, t)-sin(x + t);
% tc.u_xx = @(x, t)-cos(x + t);

% domain of computation
tc.dm = [0, 2 * pi];

end
