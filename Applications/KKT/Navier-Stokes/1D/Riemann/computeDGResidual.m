% Compute the DG residual
function res = computeDGResidual(msh, md, U, isViscous, IME, tc, bs)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(U, [bs.nb, 3 * msh.nLElems]);
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
    res = zeros(1, 3 * msh.nLElems);
else
    p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
    F = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    if isViscous
        F(bs.nep + 1 : 2 * bs.nep, :) = F(bs.nep + 1 : 2 * bs.nep, :) - tau;
        F(2 * bs.nep + 1 : end, :)    = F(2 * bs.nep + 1 : end, :)    - tau .* v + q;
    end
    
    res = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
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
alpha = max(getEigenMax(ul1(1, :), vl1, pl, tc), getEigenMax(ur1(1, :), vr1, pr, tc));

F_hat = 0.5 * (Fl + Fr - alpha .* (ur1 - ul1));
res(:, (-2 : 0)' + 3 * leLIDs1) = res(:, (-2 : 0)' + 3 * leLIDs1) + bs.phitw_face{1, 2} * reshape(F_hat, [1, 3 * nf1]);
res(:, (-2 : 0)' + 3 * reLIDs1) = res(:, (-2 : 0)' + 3 * reLIDs1) - bs.phitw_face{1, 1} * reshape(F_hat, [1, 3 * nf1]);

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
alpha = max(getEigenMax(ul(1, :), vl, pl, tc), getEigenMax(ur(1, :), vr, pr, tc));

F_hat = 0.5 * (-Fl - Fr - alpha .* (ur - ul));
res(:, (-2 : 0)' + 3 * leLIDs) = res(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 1} * F_hat';

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

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
res(:, (-2 : 0)' + 3 * leLIDs) = res(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * F_hat';

%**************************************************************************   
%                            DG residual 
%************************************************************************** 
res = res(:, (1 : 3 : end)' + (0 : 2));
res = res(:);

end