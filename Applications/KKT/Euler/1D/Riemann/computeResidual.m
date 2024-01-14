% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
u = bs.phi{1} * U;
p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
if (bs.type == 100)
    Lh = zeros(1, 3 * msh.nLElems);
else 
    F = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    Lh = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = getPressure(ur(1, :), ur(2, :), ur(3, :), tc);
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(ur(1, :), ur(2, :), ur(3, :), pr);

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));
waveVel = max(alpha);

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
Lh(:, (-2 : 0)' + 3 * leLIDs) = Lh(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * reshape(F_hat, [1, 3 * nf]);
Lh(:, (-2 : 0)' + 3 * reLIDs) = Lh(:, (-2 : 0)' + 3 * reLIDs) - bs.phitw_face{1, 1} * reshape(F_hat, [1, 3 * nf]);

%**************************************************************************   
%               boundary face contributions for DG residual
%**************************************************************************  
% left boundary face contributions   
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoL; tc.mL; tc.EL];
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pL;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(tc.rhoL, tc.mL, tc.EL, tc.pL);

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));
waveVel = max(waveVel, alpha);

F_hat = 0.5 * (-Fl - Fr - alpha .* (ur - ul));
Lh(:, (-2 : 0)' + 3 * leLIDs) = Lh(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 1} * F_hat';

% right boundary face contributions   
faceIDs = md.bndLFaces{2, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoR; tc.mR; tc.ER];
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pR;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(tc.rhoR, tc.mR, tc.ER, tc.pR);

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));
waveVel = max(waveVel, alpha);

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
Lh(:, (-2 : 0)' + 3 * leLIDs) = Lh(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * F_hat';

%**************************************************************************   
%                          residual vector
%************************************************************************** 
% Compute DG residual for current Runge-Kutta stage
LStage(:, stage.Ls) = reshape(Lh(:, (1 : 3 : end)' + (0 : 2)), [3 * KKT.NDOF, 1]);

% DG discretization contribution for residual vector
res = zeros(3 * KKT.NDOF + 2 * KKT.LDOF, 1);
if (bt.form == 1) 
    res(1 : 3 * KKT.NDOF) = M * (KStage(:, 1) - U0) + dt * (LStage(:, 1 : stage.s) * bt.A(stage.s, 1 : stage.s)');
elseif (bt.form == 2)  
    res(1 : 3 * KKT.NDOF) = M * KStage(:, stage.s) + LStage(:, 1);
end

% Add positivity constraints contribution
if (KKT.positivity)
    if (bt.form == 1) 
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;       
    end    
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 3
        dp = getPressureDeriv(v, i, tc);
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - cdt * reshape(bs.phi{1}' * (dp .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end

    g = KKT.Umin - [u(:, 1 : 3 : end), p];
    res(3 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end


