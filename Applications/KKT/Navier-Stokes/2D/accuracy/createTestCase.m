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
tc.rho0 = @(x, y)ones(size(x));
tc.u0   = @(x, y)zeros(size(x));
tc.v0   = @(x, y)zeros(size(x));
tc.p0   = @(x, y)1.e-4 + sin(x + y).^8;
tc.m10  = @(x, y)tc.rho0(x, y) .* tc.u0(x, y);
tc.m20  = @(x, y)tc.rho0(x, y) .* tc.v0(x, y);
tc.E0   = @(x, y)tc.p0(x, y) / (tc.gamma - 1) + 0.5 * tc.rho0(x, y) .* (tc.u0(x, y).^2 + tc.v0(x, y).^2);

% exact solutions and their derivatives
tc.rho = @(x, y, t)ones(size(x));
tc.u   = @(x, y, t)zeros(size(x));
tc.v   = @(x, y, t)zeros(size(x));
tc.p   = @(x, y, t)1.e-4 + sin(x + y + t).^8;

tc.rho_t  = @(x, y, t)zeros(size(x));
tc.rho_x  = @(x, y, t)zeros(size(x));
tc.rho_y  = @(x, y, t)zeros(size(x));
tc.rho_xx = @(x, y, t)zeros(size(x));
tc.rho_yy = @(x, y, t)zeros(size(x));

tc.u_t  = @(x, y, t)zeros(size(x));
tc.u_x  = @(x, y, t)zeros(size(x));
tc.u_y  = @(x, y, t)zeros(size(x));
tc.u_xx = @(x, y, t)zeros(size(x));
tc.u_xy = @(x, y, t)zeros(size(x));
tc.u_yy = @(x, y, t)zeros(size(x));

tc.v_t  = @(x, y, t)zeros(size(x));
tc.v_x  = @(x, y, t)zeros(size(x));
tc.v_y  = @(x, y, t)zeros(size(x));
tc.v_xx = @(x, y, t)zeros(size(x));
tc.v_xy = @(x, y, t)zeros(size(x));
tc.v_yy = @(x, y, t)zeros(size(x));

tc.p_t  = @(x, y, t)8 * sin(x + y + t).^7 .* cos(x + y + t);
tc.p_x  = @(x, y, t)8 * sin(x + y + t).^7 .* cos(x + y + t);
tc.p_y  = @(x, y, t)8 * sin(x + y + t).^7 .* cos(x + y + t);
tc.p_xx = @(x, y, t)8 * (7 * sin(x + y + t).^6 .* cos(x + y + t).^2 - sin(x + y + t).^8);
tc.p_yy = @(x, y, t)8 * (7 * sin(x + y + t).^6 .* cos(x + y + t).^2 - sin(x + y + t).^8);

% domain of computation
tc.dm = [0, 2 * pi, 0, 2 * pi];

end


