function Lamb = checkCompatibility(msh, U, Lamb, bs, KKT)

% Compute density at all quadrature points 
rho = bs.phi{1} * reshape(U(1 : KKT.NDOF), [bs.nb, msh.nLElems]);

% Define positivity constraints at each limiter point
g = KKT.Umin - rho(:);

% Check if the Lagrange multiplier is positive or zero
Lamb(Lamb < 0) = 0;

% If constraint is satisfied then the Lagrange multiplier must be zero 
% otherwise it must be positive and must be large than or equal to Lamb0
Lamb = (g > 0) .* max(Lamb, KKT.Lamb0);

end