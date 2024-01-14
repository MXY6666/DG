function tc = createTestCase(gamma)

% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;
        
% initial solution
tc.rho0 = @(x)1.e-7 + sin(x).^8;
tc.u0   = @(x)ones(size(x));
tc.p0   = @(x)ones(size(x));
tc.m0   = @(x)tc.rho0(x);
tc.E0   = @(x)1 / (tc.gamma - 1) + 0.5 * tc.rho0(x);

% exact solution
tc.rho = @(x, t)1.e-7 + sin(x - t).^8;     
tc.u   = @(x, t)ones(size(x));
tc.p   = @(x, t)ones(size(x));

% domain of computation
tc.dm = [0, 2 * pi];

end