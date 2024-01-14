function Lamb = checkCompatibility(msh, U0, KStage, Lamb, dt, stage, A, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);
W = A * U(1 : end / 3);

% Compute physical quantities at all quadrature points 
u = bs.phi{1} * reshape(U, [bs.nb, 3 * msh.nLElems]);
w = bs.phi{1} * reshape(W, [bs.nb, msh.nLElems]);
v = u(:, msh.nLElems + 1 : 2 * msh.nLElems) ./ u(:, 1 : msh.nLElems);
T = getTemperature(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, tc);
p = getPressure(u(:, 1 : msh.nLElems), T);

% Define positivity constraints at each limiter point
g = KKT.Umin - [u(:, 1 : msh.nLElems), p];

% Check if the Lagrange multiplier is positive or zero
Lamb(Lamb < 0) = 0;

% If constraint is satisfied then the Lagrange multiplier must be zero 
% otherwise it must be positive and must be large than or equal to Lamb0
Lamb = (g(:) > 0) .* max(Lamb, KKT.Lamb0);

end