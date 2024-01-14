function tc = createTestCase(gamma)

% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;
        
% initial solution
tc.rho_fs = 1;
tc.u_fs   = 0;
tc.v_fs   = 0;
tc.E_fs   = 1.e-9;
tc.m1_fs  = tc.rho_fs * tc.u_fs;
tc.m2_fs  = tc.rho_fs * tc.v_fs;
tc.Em     = 0.244816;

% domain of computation
tc.dm = [0, 1.1, 0, 1.1];

end