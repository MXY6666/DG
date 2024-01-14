function tc = createTestCase(gamma, Pr, Re)

% ratio of specific heats
tc.gamma = gamma;

% Prandtl number
tc.Pr = Pr;

% Reynolds number
tc.Re = Re;

% dynamic viscous coefficient
tc.mu = 1 / tc.Re;

% initial solutions
tc.rho0 = @(x)ones(size(x));  
tc.u0   = @(x)zeros(size(x));
tc.p0   = @(x)1.e-10 + sin(x).^8;
tc.m0   = @(x)tc.rho0(x) .* tc.u0(x);
tc.E0   = @(x)tc.p0(x) / (tc.gamma - 1) + 0.5 * tc.rho0(x).* tc.u0(x).^2;

% exact solutions and their derivatives
tc.rho = @(x, t)ones(size(x));  
tc.u   = @(x, t)zeros(size(x));
tc.p   = @(x, t)1.e-10 + sin(x + t).^8;

tc.rho_t  = @(x, t)zeros(size(x));
tc.rho_x  = @(x, t)zeros(size(x));
tc.rho_xx = @(x, t)zeros(size(x));

tc.u_t  = @(x, t)zeros(size(x));
tc.u_x  = @(x, t)zeros(size(x));
tc.u_xx = @(x, t)zeros(size(x));

tc.p_t  = @(x, t)8 * sin(x + t).^7 .* cos(x + t);
tc.p_x  = @(x, t)8 * sin(x + t).^7 .* cos(x + t);
tc.p_xx = @(x, t)8 * (7 * sin(x + t).^6 .* cos(x + t).^2 - sin(x + t).^8);

% domain of computation
tc.dm = [0, 2 * pi];

end


