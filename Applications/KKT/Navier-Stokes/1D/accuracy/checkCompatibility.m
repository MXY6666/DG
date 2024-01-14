function Lamb = checkCompatibility(msh, U0, KStage, Lamb, dt, stage, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);

% Compute physical quantities at all quadrature points 
u = bs.phi{1} * reshape(U, [bs.nb, 3 * msh.nLElems]);
p = getPressure(u(:, 1 : msh.nLElems), u(:, msh.nLElems + 1 : 2 * msh.nLElems), u(:, 2 * msh.nLElems + 1 : end), tc);

% Define positivity constraints at each limiter point
g = KKT.Umin - [u(:, 1 : msh.nLElems), p];

% If constraint is satisfied then the Lagrange multiplier must be zero 
% otherwise it must be positive and must be large than or equal to Lamb0
Lamb = (g(:) > 0) .* max(Lamb, KKT.Lamb0);

end