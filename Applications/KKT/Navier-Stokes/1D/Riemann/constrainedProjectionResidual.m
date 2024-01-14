function res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad1, bs, KKT)

U0 = reshape(U0, [bs.nb, 3 * msh.nLElems]);
U  = reshape(U,  [bs.nb, 3 * msh.nLElems]);
% L2 projection contribution
res = zeros(3 * KKT.NDOF + 2 * KKT.LDOF, 1);
uh  = bs.phi{1} * U;
if (flag == 1)
    u0 = zeros(quad1.np, 3 * msh.nLElems);
    if (tc.pro == 8)
        u0(:, 1 : end / 10)                    = tc.rhoL;
        u0(:, end / 10 + 1 : end / 3)          = tc.rhoR;
        u0(:, end / 3 + 1 : 13 * end / 30)     = tc.mL;
        u0(:, 13 * end / 30 + 1 : 2 * end / 3) = tc.mR;
        u0(:, 2 * end / 3 + 1 : 23 * end / 30) = tc.EL;
        u0(:, 23 * end / 30 + 1 : end)         = tc.ER;     
    else
        u0(:, 1 : end / 6)                   = tc.rhoL;
        u0(:, end / 6 + 1 : end / 3)         = tc.rhoR;
        u0(:, end / 3 + 1 : end / 2)         = tc.mL;
        u0(:, end / 2 + 1 : 2 * end / 3)     = tc.mR;
        u0(:, 2 * end / 3 + 1 : 5 * end / 6) = tc.EL;
        u0(:, 5 * end / 6 + 1 : end)         = tc.ER; 
    end  
    if (tc.pro == 10)
        u0(:, 5 * end / 6) = tc.Em;
    end
    if (tc.pro == 11)
        u0(:, 5 * end / 6) = tc.Em / msh.elemLength(:, 1);    
    end
    res(1 : 3 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 3]) .* (bs.phitw{1} * (uh - u0));
elseif (flag == 2)
    res(1 : 3 * KKT.NDOF) = repmat(msh.elemJac(:, msh.LElems), [1, 3]) .* (bs.phitw{1} * bs.phi{1} * (U - U0));
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    vh = uh(:, msh.nLElems + 1 : 2 * msh.nLElems) ./ uh(:, 1 : msh.nLElems);
    ph = getPressure(uh(:, 1 : msh.nLElems), uh(:, msh.nLElems + 1 : 2 * msh.nLElems), uh(:, 2 * msh.nLElems + 1 : end), tc);
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 3
        dp = getPressureDeriv(vh, i, tc);
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - reshape(bs.phi{1}' * (dp .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end    

    g = KKT.Umin - [uh(:, 1 : msh.nLElems), ph];
    res(3 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end