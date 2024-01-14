function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, A, tc, quad1, bs, KKT)

U0 = reshape(U0, [bs.nb, 3 * msh.nLElems]);
U  = reshape(U,  [bs.nb, 3 * msh.nLElems]);
% L2 projection contribution
res = zeros(3 * KKT.NDOF + 2 * KKT.LDOF, 1);
uh  = bs.phi{1} * U;
if (flag == 1)
    gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points;
    res(1 : 3 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 3]) .* (bs.phitw{1} * (uh - [tc.rho0(gp), tc.m0(gp), tc.E0(gp)]));
elseif (flag == 2)
    res(1 : 3 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 3]) .* (bs.phitw{1} * bs.phi{1} * (U - U0));
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    wh  = bs.phi{1} * reshape(A * U(1 : end / 3)', [bs.nb, msh.nLElems]);
    vh  = uh(:, msh.nLElems + 1 : 2 * msh.nLElems) ./ uh(:, 1 : msh.nLElems);
    Th  = getTemperature(uh(:, 1 : msh.nLElems), vh, uh(:, 2 * msh.nLElems + 1 : end), wh, tc);
    ph  = getPressure(uh(:, 1 : msh.nLElems), Th);
    dph = getPressureDeriv(uh(:, 1 : msh.nLElems), vh, uh(:, 2 * msh.nLElems + 1 : end), wh, tc);
       
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 3
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - reshape(bs.phi{1}' * (dph(:, :, i) .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end    
    tmp = bs.phi{1}' * (dph(:, :, 4) .* Lamb(:, msh.nLElems + 1 : end));
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - (A' * tmp(:));

    g = KKT.Umin - [uh(:, 1 : msh.nLElems), ph];
    res(3 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end