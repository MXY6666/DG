function tc = createTestCase(Re, We, Pr, Cv)

% Reynolds number
tc.Re = Re;

% Weber number
tc.We = We;

% Prandtl number
tc.Pr = Pr;

% specific heat at constant volume
tc.Cv = Cv;

% initial solution
tc.rho0 = @(x)ones(size(x));
tc.u0   = @(x)zeros(size(x));
tc.T0   = @(x)0.75 + 0.25 * (1.e-10 + sin(x).^8);
rhox0   = @(x)zeros(size(x));
tc.m0   = @(x)tc.rho0(x) .* tc.u0(x);
tc.E0   = @(x)tc.rho0(x) .* (8 / 3 * tc.Cv * tc.T0(x) - 3 * tc.rho0(x) + 0.5 * tc.u0(x).^2) + 0.5 / tc.We * rhox0(x).^2;

% exact solutions and their derivatives
tc.rho = @(x, t)ones(size(x));
tc.u   = @(x, t)zeros(size(x));
tc.T   = @(x, t)0.75 + 0.25 * (1.e-10 + sin(x + t).^8);
tc.p   = @(x, t)1.e-10 + sin(x + t).^8;

tc.rho_t   = @(x, t)zeros(size(x));
tc.rho_x   = @(x, t)zeros(size(x)); 
tc.rho_xx  = @(x, t)zeros(size(x));
tc.rho_xt  = @(x, t)zeros(size(x));
tc.rho_xxx = @(x, t)zeros(size(x));

tc.u_t  = @(x, t)zeros(size(x));
tc.u_x  = @(x, t)zeros(size(x));
tc.u_xx = @(x, t)zeros(size(x));

tc.T_t  = @(x, t)2 * sin(x + t).^7 .* cos(x + t);
tc.T_x  = @(x, t)2 * sin(x + t).^7 .* cos(x + t);
tc.T_xx = @(x, t)2 * (7 * sin(x + t).^6 .* cos(x + t).^2 - sin(x + t).^8);

% domain of computation
tc.dm = [0, 2 * pi];

end


