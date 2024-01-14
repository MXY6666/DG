function [U0, Lamb] = computeInitialSolution(msh, md, tc, quad, bs, slo, KKT, newton, mtol)

% initial solution
U0 = zeros(bs.nb, 3 * msh.nLElems);
if (tc.pro == 8)
    U0(1, 1 : end / 10)                    = tc.rhoL;
    U0(1, end / 10 + 1 : end / 3)          = tc.rhoR;
    U0(1, end / 3 + 1 : 13 * end / 30)     = tc.mL;
    U0(1, 13 * end / 30 + 1 : 2 * end / 3) = tc.mR;
    U0(1, 2 * end / 3 + 1 : 23 * end / 30) = tc.EL;
    U0(1, 23 * end / 30 + 1 : end)         = tc.ER;  
else
    U0(1, 1 : end / 6)                   = tc.rhoL;
    U0(1, end / 6 + 1 : end / 3)         = tc.rhoR;
    U0(1, end / 3 + 1 : end / 2)         = tc.mL;
    U0(1, end / 2 + 1 : 2 * end / 3)     = tc.mR;
    U0(1, 2 * end / 3 + 1 : 5 * end / 6) = tc.EL;
    U0(1, 5 * end / 6 + 1 : end)         = tc.ER;  
end
if (tc.pro == 10)
    U0(1, 5 * end / 6) = tc.Em; 
end
if (tc.pro == 11)
    U0(1, 5 * end / 6) = tc.Em / msh.elemLength(:, 1); 
end
U0 = U0(:);

% Apply slope limiter to suppress oscillation
U0 = applySlopeLimiter(msh, md, U0, tc, bs, slo);

% Apply KKT limiter to preserve positivity
if (KKT.enable)
    [U0, Lamb, fail] = constrainedProjection(msh, U0, 0, 0, 1, tc, quad, bs, KKT, newton, mtol);
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
