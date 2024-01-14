function [U0, Lamb] = computeInitialSolution(msh, tc, quad1, bs, KKT, newton, mtol)

U0 = computeInitialSolution_line(msh, {tc.rho0, tc.m0, tc.E0}, bs);
U0 = U0(:);
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, tc, quad1, bs, KKT, newton, mtol);
    if (fail)
       error('Constrained projection for initial solution was not feasible')
    end
else
    Lamb = [];
end
u = bs.phi{1} * reshape(U0, [bs.nb, 3 * msh.nLElems]);
p = getPressure(u(:, 1 : msh.nLElems), u(:, msh.nLElems + 1 : 2 * msh.nLElems), u(:, 2 * msh.nLElems + 1 : end), tc);
rhomin = min(u(:, 1 : msh.nLElems), [], 'all');
pmin   = min(p, [], 'all');
fprintf('\n\nExtrema at time          %+15.6e\n', 0);
fprintf('Minimum density          %+15.6e\n', rhomin);
fprintf('Minimum pressure         %+15.6e\n', pmin);

end