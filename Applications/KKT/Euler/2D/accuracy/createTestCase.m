function tc = createTestCase(gamma)

% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;
        
% initial solution
tc.rho0 = @(x, y)2 + sin(2 * pi * (x + y));
tc.u0   = @(x, y)ones(size(x));
tc.v0   = @(x, y)ones(size(x));
tc.p0   = @(x, y)ones(size(x));
tc.m10  = @(x, y)tc.rho0(x, y);
tc.m20  = @(x, y)tc.rho0(x, y);
tc.E0   = @(x, y)1 ./ (tc.gamma - 1) + tc.rho0(x, y);

% exact solution
tc.rho = @(x, y, t)2 + sin(2 * pi * (x + y - 2 * t));
tc.u   = @(x, y, t)ones(size(x));
tc.v   = @(x, y, t)ones(size(x));
tc.p   = @(x, y, t)ones(size(x));

% domain of computation
tc.dm = [0, 1, 0, 1];

end