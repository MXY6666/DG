function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad, bs, KKT)

U0 = reshape(U0, [bs.nb, msh.nLElems]);
U  = reshape(U,  [bs.nb, msh.nLElems]);
% L2 projection contribution
res = zeros(KKT.NDOF + KKT.LDOF, 1);
uh = bs.phi{1} * U;
if (flag == 1)
    gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad.points;
    res(1 : KKT.NDOF) = msh.elemJac(:, msh.LElems) .* (bs.phitw{1} * (uh - tc.u0(gp)));
elseif (flag == 2)
    res(1 : KKT.NDOF) = msh.elemJac(:, msh.LElems) .* (bs.phitw{1} * bs.phi{1} * (U - U0));
end

% Add contribution of constraint equations
if (KKT.enable)
    inc = KKT.positivity + KKT.maximum;
    Lamb = reshape(Lamb, [bs.nep, inc * msh.nLElems]);
    
    % Add positivity constraints to DG discretization
    if (KKT.positivity)
        res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb(:, 1 : inc : end), [KKT.NDOF, 1]);
    end
    
    % Add Maximum principle constraints to DG discretization
    if (KKT.maximum)
        res(1 : KKT.NDOF) = res(1 : KKT.NDOF) + reshape(bs.phi{1}' * Lamb(:, KKT.positivity + 1 : inc : end), [KKT.NDOF, 1]); 
    end
    
    % Add contributions from inequality constraints
    g = zeros(bs.nep, inc * msh.nLElems);    
    if (KKT.positivity)
        g(:, 1 : inc : end) = KKT.Umin - uh;
    end
    if (KKT.maximum)
        g(:, KKT.positivity + 1 : inc : end) = uh - KKT.Umax;
    end 
    
    res(KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end