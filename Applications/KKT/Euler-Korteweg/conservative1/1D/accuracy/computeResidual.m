function res = computeResidual(msh, U, Lamb, AU0, res1, k, m, dts, bs, SDC, KKT)

res = zeros(2 * KKT.NDOF + KKT.LDOF, 1);
% DG discretization contribution
if (k == 1)
    res(1 : 2 * KKT.NDOF) = res1 + U - dts(m - 1) * (AU0 * U);        
else
    res(1 : 2 * KKT.NDOF) = res1 + U - SDC.theta * dts(m - 1) * (AU0 * U);
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, msh.nLElems]);
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb, [KKT.NDOF, 1]);

    rho = bs.phi{1} * reshape(U(1 : KKT.NDOF), [bs.nb, msh.nLElems]);
    g = KKT.Umin - rho(:);
    res(2 * KKT.NDOF + 1 : end) = min(-g, Lamb(:));
end

end