function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad2, bs, KKT)

U0 = reshape(U0, [bs.nb, msh.nLElems]);
U  = reshape(U,  [bs.nb, msh.nLElems]);
% L2 projection contribution
res = zeros(KKT.NDOF + KKT.LDOF, 1);
uh  = bs.phi{1} * U;
if (flag == 1)
    ct  = msh.elemCenter(:, msh.LElems);
    h   = msh.elemLength(:, msh.LElems);
    gpx = ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1);
    gpy = ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2);
    res(1 : KKT.NDOF) = msh.elemJac(:, msh.LElems) .* (bs.phitw{1} * (uh - tc.u0(gpx, gpy)));
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