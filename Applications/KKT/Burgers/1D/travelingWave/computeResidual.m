function [res, LStage] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);
U = reshape(U, [bs.nb, msh.nLElems]);

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
Lh = -bs.phitw{2} * f(bs.phi{1} * U);

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************
faceIDs = md.intLFaces{1, 3};
leLIDs = msh.faceElems(1, faceIDs);
reLIDs = msh.faceElems(2, faceIDs);

ul = bs.phi_face{1, 2} * U(:, leLIDs); 
ur = bs.phi_face{1, 1} * U(:, reLIDs);
fl = f(ul); 
fr = f(ur);

alpha = max(abs(ul), abs(ur)); 
f_hat = 0.5 * (fl+ fr + alpha .* (ul - ur));
Lh(:, leLIDs) = Lh(:, leLIDs) + bs.phitw_face{1, 2} * f_hat;
Lh(:, reLIDs) = Lh(:, reLIDs) - bs.phitw_face{1, 1} * f_hat;

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
