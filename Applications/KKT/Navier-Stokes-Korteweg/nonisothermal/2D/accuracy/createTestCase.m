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
tc.rho0 = @(x, y)ones(size(x));
tc.u0   = @(x, y)zeros(size(x));
tc.v0   = @(x, y)zeros(size(x));
tc.T0   = @(x, y)0.75 + 0.25 * (1.e-10 + sin(x + y).^8);
rhox0   = @(x, y)zeros(size(x));
rhoy0   = @(x, y)zeros(size(x));
tc.m10  = @(x, y)tc.rho0(x, y) .* tc.u0(x, y);
tc.m20  = @(x, y)tc.rho0(x, y) .* tc.v0(x, y);
tc.E0   = @(x, y)tc.rho0(x, y) .* (8 / 3 * tc.Cv * tc.T0(x, y) - 3 * tc.rho0(x, y) + 0.5 * (tc.u0(x, y).^2 + tc.v0(x, y).^2)) + 0.5 / tc.We * (rhox0(x, y).^2 + rhoy0(x, y).^2);

% exact solutions and their derivatives
tc.rho = @(x, y, t)ones(size(x));
tc.u   = @(x, y, t)zeros(size(x));
tc.v   = @(x, y, t)zeros(size(x));
tc.T   = @(x, y, t)0.75 + 0.25 * (1.e-10 + sin(x + y + t).^8);
tc.p   = @(x, y, t)1.e-10 + sin(x + y + t).^8;

tc.rho_x   = @(x, y, t)zeros(size(x));
tc.rho_y   = @(x, y, t)zeros(size(x));
tc.rho_t   = @(x, y, t)zeros(size(x));
tc.rho_xx  = @(x, y, t)zeros(size(x));
tc.rho_yy  = @(x, y, t)zeros(size(x));
tc.rho_xy  = @(x, y, t)zeros(size(x));
tc.rho_xt  = @(x, y, t)zeros(size(x));
tc.rho_yt  = @(x, y, t)zeros(size(x));
tc.rho_xxx = @(x, y, t)zeros(size(x));
tc.rho_yyy = @(x, y, t)zeros(size(x));
tc.rho_xxy = @(x, y, t)zeros(size(x));
tc.rho_xyy = @(x, y, t)zeros(size(x));

tc.u_x  = @(x, y, t)zeros(size(x));
tc.u_y  = @(x, y, t)zeros(size(x));
tc.u_t  = @(x, y, t)zeros(size(x));
tc.u_xx = @(x, y, t)zeros(size(x));
tc.u_xy = @(x, y, t)zeros(size(x));
tc.u_yy = @(x, y, t)zeros(size(x));

tc.v_x  = @(x, y, t)zeros(size(x));
tc.v_y  = @(x, y, t)zeros(size(x));
tc.v_t  = @(x, y, t)zeros(size(x));
tc.v_xx = @(x, y, t)zeros(size(x));
tc.v_xy = @(x, y, t)zeros(size(x));
tc.v_yy = @(x, y, t)zeros(size(x));

tc.T_x  = @(x, y, t)2 * sin(x + y + t).^7 .* cos(x + y + t);
tc.T_y  = @(x, y, t)2 * sin(x + y + t).^7 .* cos(x + y + t);
tc.T_t  = @(x, y, t)2 * sin(x + y + t).^7 .* cos(x + y + t);
tc.T_xx = @(x, y, t)2 * (7 * sin(x + y + t).^6 .* cos(x + y + t).^2 - sin(x + y + t).^8);
tc.T_yy = @(x, y, t)2 * (7 * sin(x + y + t).^6 .* cos(x + y + t).^2 - sin(x + y + t).^8);

% domain of computation
tc.dm = [0, 2 * pi, 0, 2 * pi];

end

