function Lamb = checkCompatibility(msh, U0, KStage, Lamb, dt, stage, bs, bt, KKT)

inc = KKT.positivity + KKT.maximum;
% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);

% Compute solution at all limiter points 
u = bs.phi{1} * reshape(U, [bs.nb, msh.nLElems]);

% Compute positivity constraints at each limiter point
g = zeros(bs.nep, inc * msh.nLElems);
if (KKT.positivity)
    g(:, 1 : inc : end) = KKT.Umin - u;
end
if (KKT.maximum)
    g(:, KKT.positivity + 1 : inc : end) = u - KKT.Umax;
end  
g = g(:);

% If constraint is satisfied then the Lagrange multiplier must be zero 
% otherwise it must be positive and must be large than or equal to Lamb0
Lamb = (g > 0) .* max(Lamb, KKT.Lamb0);

end