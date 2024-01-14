function Lamb = checkCompatibility(msh, U0, KStage, Lamb, dt, stage, Ax, Ay, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U  = getStageSolution(U0, KStage, dt, stage, bt);
W1 = Ax * U(1 : end / 4);
W2 = Ay * U(1 : end / 4);

% Compute physical quantities at all quadrature points 
u  = bs.phi{1} * reshape(U, [bs.nb, 4 * msh.nLElems]);
w1 = bs.phi{1} * reshape(W1, [bs.nb, msh.nLElems]);
w2 = bs.phi{1} * reshape(W2, [bs.nb, msh.nLElems]);
v  = u(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(u(:, 1 : msh.nLElems), [1, 2]);
T  = getTemperature(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, tc);
p  = getPressure(u(:, 1 : msh.nLElems), T);

% Define positivity constraints at each limiter point
g = KKT.Umin - [u(:, 1 : msh.nLElems), p];

% Check if the Lagrange multiplier is positive or zero
Lamb(Lamb < 0) = 0;

% If constraint is satisfied then the Lagrange multiplier must be zero 
% otherwise it must be positive and must be large than or equal to Lamb0
Lamb = (g(:) > 0) .* max(Lamb, KKT.Lamb0);

end