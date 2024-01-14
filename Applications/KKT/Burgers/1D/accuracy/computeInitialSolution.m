function [U0, Lamb] = computeInitialSolution(msh, tc, quad, bs, KKT, newton)

U0 = computeInitialSolution_line(msh, {tc.u0}, bs);
U0 = U0(:);
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, tc, quad, bs, KKT, newton);
    if (fail)
       error('Constrained projection for initial solution was not feasible')
    end
else
    Lamb = [];
end
u = bs.phi{1} * reshape(U0, [bs.nb, msh.nLElems]);
umin = min(u, [], 'all');
umax = max(u, [], 'all');
fprintf('\n\nExtrema at time          %+15.6e\n', 0);
fprintf('Minimum uh               %+15.6e\n', umin);
fprintf('Maximum uh               %+15.6e\n', umax);


end