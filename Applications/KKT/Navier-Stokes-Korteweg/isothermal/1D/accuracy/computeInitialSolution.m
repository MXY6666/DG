function [U0, Lamb] = computeInitialSolution(msh, tc, quad1, bs, KKT, newton)

U0 = computeInitialSolution_line(msh, {tc.rho0, tc.m0}, bs);
U0 = U0(:);
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, tc, quad1, bs, KKT, newton);
    if (fail)
       error('Constrained projection for initial solution was not feasible')
    end
else
    Lamb = [];
end
rho = bs.phi{1} * reshape(U0(1 : KKT.NDOF), [bs.nb, msh.nLElems]);
rhomin = min(rho, [], 'all');
fprintf('\n\nExtrema at time          %+15.6e\n', 0);
fprintf('Minimum density          %+15.6e\n', rhomin);

end