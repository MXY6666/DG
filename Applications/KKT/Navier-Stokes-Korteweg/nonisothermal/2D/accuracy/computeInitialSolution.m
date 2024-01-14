function [U0, W10, W20, Lamb] = computeInitialSolution(msh, Ax, Ay, tc, quad2, bs, KKT, newton, mtol)

U0 = computeInitialSolution_rect(msh, {tc.rho0, tc.m10, tc.m20, tc.E0}, bs);
U0 = U0(:);
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, Ax, Ay, tc, quad2, bs, KKT, newton, mtol);
    if (fail)
       error('Constrained projection for initial solution was not feasible')
    end
else
    Lamb = [];
end
W10  = Ax * U0(1 : end / 4); 
W20  = Ay * U0(1 : end / 4);

u  = bs.phi{1} * reshape(U0, [bs.nb, 4 * msh.nLElems]);
w1 = bs.phi{1} * reshape(W10, [bs.nb, msh.nLElems]);
w2 = bs.phi{1} * reshape(W20, [bs.nb, msh.nLElems]);
v  = u(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(u(:, 1 : msh.nLElems), [1, 2]);
T  = getTemperature(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, tc);
p  = getPressure(u(:, 1 : msh.nLElems), T);
rhomin = min(u(:, 1 : msh.nLElems), [], 'all');
pmin   = min(p, [], 'all');
fprintf('\n\nExtrema at time          %+15.6e\n', 0);
fprintf('Minimum density          %+15.6e\n', rhomin);
fprintf('Minimum pressure         %+15.6e\n', pmin);

end