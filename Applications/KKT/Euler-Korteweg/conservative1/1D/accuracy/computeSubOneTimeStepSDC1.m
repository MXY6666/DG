% Compute one time step in a time subinterval 
function [U, AU, iL, L] = computeSubOneTimeStepSDC1(msh, U0, AU0, iL0, L0, k, m, ts, dts, flag, A11, A12, Acap, IME, tc, quad1, bs, SDC, newton)
    
NDOF = msh.nLElems * bs.nb;
if (k == 1)
    A = speye(2 * NDOF) - dts(m - 1) * AU0;
    b = U0 + dts(m - 1) * iL0;    
else
    A = speye(2 * msh.nLElems * bs.nb) - SDC.theta * dts(m - 1) * AU0;
    b = U0 + SDC.theta * dts(m - 1) * (iL0 - L0(:, m)) + sum(dts) * sum(SDC.weights(m - 1, :) .* L0, 2); 
end

if (newton.linearSolver.type == 1)
    U = A \ b;  
elseif (newton.linearSolver.type == 3)
    [PL, PU] = ilu(A, newton.linearSolver.ilu);
    U = gmres(A, b, newton.linearSolver.restart, newton.linearSolver.tolres, newton.linearSolver.maxit, PL, PU, U0);     
end 

if (flag == 0)
    AU = sparse([], [], [], 2 * NDOF, 2 * NDOF); 
    iL = zeros(2 * NDOF, 1); 
    L  = zeros(2 * NDOF, 1);
else
    [AU, iL] = computeTimeDerivative(msh, U, ts(m - 1) + dts(m - 1), A11, A12, Acap, IME, tc, quad1, bs);
    L = AU * U + iL; 
end

end