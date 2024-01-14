function [U0, Lamb] = computeInitialSolution(msh, md, tc, quad2, bs, slo, KKT, newton, mtol)

% initial solution
U0 = zeros(bs.nb, 4 * msh.nLElems);
U0(1, 1 : end / 4)               = tc.rho_fs;
U0(1, end / 4 + 1 : end / 2)     = tc.m1_fs;
U0(1, end / 2 + 1 : 3 * end / 4) = tc.m2_fs;
U0(1, 3 * end / 4 + 1 : end)     = tc.E_fs;
U0(1, 3 * end / 4 + 1)           = tc.Em / msh.elemSize(:, 1);
U0 = U0(:);

% Apply slope limiter to suppress oscillation
U0 = applySlopeLimiter(msh, md, U0, tc, bs, slo);

% Apply KKT limiter to preserve positivity
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, tc, quad2, bs, KKT, newton, mtol);
    if (fail)
       error('Constrained projection for initial solution was not feasible')
    end
else
    Lamb = [];
end
u = bs.phi{1} * reshape(U0, [bs.nb, 4 * msh.nLElems]);
p = getPressure(u(:, 1 : msh.nLElems), u(:, msh.nLElems + 1 : 2 * msh.nLElems), u(:, 2 * msh.nLElems + 1 : 3 * msh.nLElems), u(:, 3 * msh.nLElems + 1 : end), tc);
rhomin = min(u(:, 1 : msh.nLElems), [], 'all');
pmin   = min(p, [], 'all');
fprintf('\n\nExtrema at time          %+15.6e\n', 0);
fprintf('Minimum density          %+15.6e\n', rhomin);
fprintf('Minimum pressure         %+15.6e\n', pmin);

end