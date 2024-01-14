function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad1, bs, KKT)

U0 = reshape(U0, [bs.nb, 2 * msh.nLElems]);
U  = reshape(U,  [bs.nb, 2 * msh.nLElems]);
% L2 projection contribution
res = zeros(2 * KKT.NDOF + KKT.LDOF, 1);
uh  = bs.phi{1} * U;
if (flag == 1)
    gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points;
    res(1 : 2 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 2]) .* (bs.phitw{1} * (uh - [tc.rho0(gp), tc.u0(gp)]));
elseif (flag == 2)
    res(1 : 2 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 2]) .* (bs.phitw{1} * bs.phi{1} * (U - U0));
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, msh.nLElems]);
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb, [KKT.NDOF, 1]);

    g = KKT.Umin - uh(:, 1 : msh.nLElems);
    res(2 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end