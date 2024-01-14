function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad2, bs, KKT)

U0 = reshape(U0, [bs.nb, 4 * msh.nLElems]);
U  = reshape(U,  [bs.nb, 4 * msh.nLElems]);
% L2 projection contribution
res = zeros(4 * KKT.NDOF + 2 * KKT.LDOF, 1);
uh  = bs.phi{1} * U;
if (flag == 1)
    u0 = zeros(quad2.np, 4 * msh.nLElems);
    u0(:, 1 : end / 4)               = tc.rho_fs;
    u0(:, end / 4 + 1 : end / 2)     = tc.m1_fs;
    u0(:, end / 2 + 1 : 3 * end / 4) = tc.m2_fs;
    u0(:, 3 * end / 4 + 1 : end)     = tc.E_fs;
    u0(:, 3 * end / 4 + 1)           = tc.Em / msh.elemSize(:, 1);
    res(1 : 4 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 4]) .* (bs.phitw{1} * (uh - u0));
elseif (flag == 2)
    res(1 : 4 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 4]) .* (bs.phitw{1} * bs.phi{1} * (U - U0));
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    vh = uh(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(uh(:, 1 : msh.nLElems), [1, 2]);
    ph = getPressure(uh(:, 1 : msh.nLElems), uh(:, msh.nLElems + 1 : 2 * msh.nLElems), uh(:, 2 * msh.nLElems + 1 : 3 * msh.nLElems), uh(:, 3 * msh.nLElems + 1 : end), tc);
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 4
        dp = getPressureDeriv(vh(:, 1 : end / 2), vh(:, end / 2 + 1 : end), i, tc);
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - reshape(bs.phi{1}' * (dp .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end    

    g = KKT.Umin - [uh(:, 1 : msh.nLElems), ph];
    res(4 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end