% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, isViscous, M, IME, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

% Evaluate physical quantities 
u = bs.phi{1} * U; 
v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);

faceIDs1 = md.intLFaces{1, 3};
leLIDs1  = msh.faceElems(1, faceIDs1);
reLIDs1  = msh.faceElems(2, faceIDs1);
nf1      = length(faceIDs1);

ul1 = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs1), [3, nf1]);
ur1 = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs1), [3, nf1]);
vl1 = ul1(2, :) ./ ul1(1, :);
vr1 = ur1(2, :) ./ ur1(1, :);

%**************************************************************************
%                     Compute viscous stress tensor
%************************************************************************** 
if isViscous
    % element contributions
    Z = -bs.phitw{2} * v;

    % internal face contributions
    v_hat = 0.5 * (vl1 + vr1);
    Z(:, leLIDs1) = Z(:, leLIDs1) + bs.phitw_face{1, 2} * v_hat;
    Z(:, reLIDs1) = Z(:, reLIDs1) - bs.phitw_face{1, 1} * v_hat;
    
    % left boundary face contributions  
    faceIDs = md.bndLFaces{1, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    
    Z(:, leLIDs) = Z(:, leLIDs) - bs.phitw_face{1, 1} * tc.uL;
    
    % right boundary face contributions 
    faceIDs = md.bndLFaces{2, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    
    Z(:, leLIDs) = Z(:, leLIDs) + bs.phitw_face{1, 2} * tc.uR;
    
    % Take care of mass mtrix
    Z = (IME * Z) ./ msh.elemJac(:, msh.LElems);

    % viscous stress tensor
    TAU = 4 / 3 * tc.mu * Z;

    % Evaluate viscous stress tensor
    tau  = bs.phi{1} * TAU;
    taul = bs.phi_face{1, 2} * TAU(:, leLIDs1);
    taur = bs.phi_face{1, 1} * TAU(:, reLIDs1);
end

%**************************************************************************
%                           Compute heat flux
%************************************************************************** 
if isViscous
    % element contributions
    Q = -bs.phitw{2} * getScaledTemperature(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), tc);

    % internal face contributions
    sT_hat = 0.5 * (getScaledTemperature(ul1(1, :), vl1, ul1(3, :), tc) + getScaledTemperature(ur1(1, :), vr1, ur1(3, :), tc));
    Q(:, leLIDs1) = Q(:, leLIDs1) + bs.phitw_face{1, 2} * sT_hat;
    Q(:, reLIDs1) = Q(:, reLIDs1) - bs.phitw_face{1, 1} * sT_hat;
    
    % left boundary face contributions  
    faceIDs = md.bndLFaces{1, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    
    Q(:, leLIDs) = Q(:, leLIDs) - bs.phitw_face{1, 1} * (-tc.kTL);
    
    % right boundary face contributions 
    faceIDs = md.bndLFaces{2, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    
    Q(:, leLIDs) = Q(:, leLIDs) + bs.phitw_face{1, 2} * (-tc.kTR);    

    % Take care of mass mtrix
    Q = (IME * Q) ./ msh.elemJac(:, msh.LElems);

    % Evaluate heat flux
    q  = bs.phi{1} * Q;
    ql = bs.phi_face{1, 2} * Q(:, leLIDs1);
    qr = bs.phi_face{1, 1} * Q(:, reLIDs1);
end

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
if (bs.type == 100)
    Lh = zeros(1, 3 * msh.nLElems);
else
    p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
    F = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    if isViscous
        F(bs.nep + 1 : 2 * bs.nep, :) = F(bs.nep + 1 : 2 * bs.nep, :) - tau;
        F(2 * bs.nep + 1 : end, :)    = F(2 * bs.nep + 1 : end, :)    - tau .* v + q;
    end
    
    Lh = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
pl = getPressure(ul1(1, :), ul1(2, :), ul1(3, :), tc);
pr = getPressure(ur1(1, :), ur1(2, :), ur1(3, :), tc);
Fl = getF(ul1(1, :), ul1(2, :), ul1(3, :), pl); 
Fr = getF(ur1(1, :), ur1(2, :), ur1(3, :), pr);
if isViscous
    Fl(2, :) = Fl(2, :) - taul;
    Fr(2, :) = Fr(2, :) - taur;
    Fl(3, :) = Fl(3, :) - taul .* vl1 + ql;
    Fr(3, :) = Fr(3, :) - taur .* vr1 + qr;
end

% dissipation coefficient 
alpha   = max(getEigenMax(ul1(1, :), vl1, pl, tc), getEigenMax(ur1(1, :), vr1, pr, tc));
waveVel = max(alpha);

F_hat = 0.5 * (Fl + Fr - alpha .* (ur1 - ul1));
Lh(:, (-2 : 0)' + 3 * leLIDs1) = Lh(:, (-2 : 0)' + 3 * leLIDs1) + bs.phitw_face{1, 2} * reshape(F_hat, [1, 3 * nf1]);
Lh(:, (-2 : 0)' + 3 * reLIDs1) = Lh(:, (-2 : 0)' + 3 * reLIDs1) - bs.phitw_face{1, 1} * reshape(F_hat, [1, 3 * nf1]);

%**************************************************************************   
%               boundary face contributions for DG residual
%**************************************************************************  
% left boundary face contributions   
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoL; tc.mL; tc.EL];
vl = ul(2, :) ./ ul(1, :);
vr = tc.uL;
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pL;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(tc.rhoL, tc.mL, tc.EL, tc.pL);
if isViscous
    taul = bs.phi_face{1, 1} * TAU(:, leLIDs);
    ql   = bs.phi_face{1, 1} * Q(:, leLIDs);
    Fl(2, :) = Fl(2, :) - taul;
    Fr(2, :) = Fr(2, :) - taul;
    Fl(3, :) = Fl(3, :) - taul .* vr + ql;
    Fr(3, :) = Fr(3, :) - taul .* vr + ql;
end

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), vl, pl, tc), getEigenMax(ur(1, :), vr, pr, tc));
waveVel = max(waveVel, alpha);

F_hat = 0.5 * (-Fl - Fr - alpha .* (ur - ul));
Lh(:, (-2 : 0)' + 3 * leLIDs) = Lh(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 1} * F_hat';

% right boundary face contributions   
faceIDs = md.bndLFaces{2, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoR; tc.mR; tc.ER];
vl = ul(2, :) ./ ul(1, :);
vr = tc.uR;
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pR;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(tc.rhoR, tc.mR, tc.ER, tc.pR);
if isViscous
    taul = bs.phi_face{1, 1} * TAU(:, leLIDs);
    ql   = bs.phi_face{1, 1} * Q(:, leLIDs);
    Fl(2, :) = Fl(2, :) - taul;
    Fr(2, :) = Fr(2, :) - taul;
    Fl(3, :) = Fl(3, :) - taul .* vr + ql;
    Fr(3, :) = Fr(3, :) - taul .* vr + ql;
end

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), vl, pl, tc), getEigenMax(ur(1, :), vr, pr, tc));
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
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    if (bt.form == 1) 
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;       
    end
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 3
        dp = getPressureDeriv(v, i, tc);
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - cdt * reshape(bs.phi{1}' * (dp .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end

    g = KKT.Umin - [u(:, 1 : 3 : end), p];
    res(3 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end
