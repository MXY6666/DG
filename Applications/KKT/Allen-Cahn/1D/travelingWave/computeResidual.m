function [res, LStage] = computeResidual(msh, md, U0, KStage, Lamb, LStage, t, dt, stage, M, A1, S2, IME, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);
U = reshape(U, [bs.nb, msh.nLElems]);

% Compute Auxiliary variable
Q = computeAuxiliaryVariable(msh, md, U, t + bt.c(stage.s) * dt, A1, IME, tc, bs);

% non-stiffness contributions for DG residual
u = bs.phi{1} * U;

Lh = bs.phitw{1} * (msh.elemJac(:, msh.LElems) .* f(u));

% stiffness contributions for DG residual 
Lh = Lh(:) + S2 * Q(:);

%**************************************************************************   
%                          residual vector
%************************************************************************** 
% Compute DG residual for current Runge-Kutta stage
LStage(:, stage.Ls) = Lh(:);

% DG discretization contribution for residual vector
res = zeros(KKT.NDOF + KKT.LDOF, 1);
if (bt.form == 1) 
    res(1 : KKT.NDOF) = M * (KStage(:, 1) - U0) + dt * (LStage(:, 1 : stage.s) * bt.A(stage.s, 1 : stage.s)');
elseif (bt.form == 2)  
    res(1 : KKT.NDOF) = M * KStage(:, stage.s) + LStage(:, 1);
end

% Add contribution of constraint equations
if (KKT.enable)
    inc = KKT.positivity + KKT.maximum;
    if (bt.form == 1) 
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;       
    end
    Lamb = reshape(Lamb, [bs.nep, inc * msh.nLElems]);
    u = bs.phi{1} * reshape(U, [bs.nb, msh.nLElems]);
      
    % Add positivity constraints to DG discretization
    if (KKT.positivity)
        res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : inc : end), [KKT.NDOF, 1]);       
    end
    
    % Add Maximum principle constraints to DG discretization
    if (KKT.maximum)
        res(1 : KKT.NDOF) = res(1 : KKT.NDOF) + cdt * reshape(bs.phi{1}' * Lamb(:, KKT.positivity + 1 : inc : end), [KKT.NDOF, 1]); 
    end

    % Add contributions from inequality constraints
    g = zeros(bs.nep, inc * msh.nLElems);    
    if (KKT.positivity)
        g(:, 1 : inc : end) = KKT.Umin - u;
    end
    if (KKT.maximum)
        g(:, KKT.positivity + 1 : inc : end) = u - KKT.Umax;
    end   
    
    res(KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end
