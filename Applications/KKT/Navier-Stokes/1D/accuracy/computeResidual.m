% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, t, dt, stage, M, IME, tc, quad1, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

% Evaluate physical quantities 
u = bs.phi{1} * U; 
v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);
p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);

faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
vl = ul(2, :) ./ ul(1, :);
vr = ur(2, :) ./ ur(1, :);
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = getPressure(ur(1, :), ur(2, :), ur(3, :), tc);

%**************************************************************************
%                     Compute viscous stress tensor
%************************************************************************** 
% element contributions
Z = -bs.phitw{2} * v;

% internal face contributions
v_hat = 0.5 * (vl + vr);
Z(:, leLIDs) = Z(:, leLIDs) + bs.phitw_face{1, 2} * v_hat;
Z(:, reLIDs) = Z(:, reLIDs) - bs.phitw_face{1, 1} * v_hat;

% Take care of mass mtrix
Z = (IME * Z) ./ msh.elemJac(:, msh.LElems);

% viscous stress tensor
TAU = 4 / 3 * tc.mu * Z;

% Evaluate viscous stress tensor
tau  = bs.phi{1} * TAU;
taul = bs.phi_face{1, 2} * TAU(:, leLIDs);
taur = bs.phi_face{1, 1} * TAU(:, reLIDs);

%**************************************************************************
%                           Compute heat flux
%************************************************************************** 
% element contributions
Q = -bs.phitw{2} * getScaledTemperature(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), tc);

% internal face contributions
sT_hat = 0.5 * (getScaledTemperature(ul(1, :), vl, ul(3, :), tc) + getScaledTemperature(ur(1, :), vr, ur(3, :), tc));
Q(:, leLIDs) = Q(:, leLIDs) + bs.phitw_face{1, 2} * sT_hat;
Q(:, reLIDs) = Q(:, reLIDs) - bs.phitw_face{1, 1} * sT_hat;

% Take care of mass mtrix
Q = (IME * Q) ./ msh.elemJac(:, msh.LElems);

% Evaluate heat flux
q  = bs.phi{1} * Q;
ql = bs.phi_face{1, 2} * Q(:, leLIDs);
qr = bs.phi_face{1, 1} * Q(:, reLIDs);

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
if (bs.type == 100)
    Lh = zeros(1, 3 * msh.nLElems);
else
    F = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    F(bs.nep + 1 : 2 * bs.nep, :) = F(bs.nep + 1 : 2 * bs.nep, :) - tau;
    F(2 * bs.nep + 1 : end, :)    = F(2 * bs.nep + 1 : end, :)    - tau .* v + q;
    
    Lh = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(ur(1, :), ur(2, :), ur(3, :), pr);
Fl(2, :) = Fl(2, :) - taul;
Fr(2, :) = Fr(2, :) - taur;
Fl(3, :) = Fl(3, :) - taul .* vl + ql;
Fr(3, :) = Fr(3, :) - taur .* vr + qr;

% dissipation coefficient 
alpha   = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));
waveVel = max(alpha);

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
Lh(:, (-2 : 0)' + 3 * leLIDs) = Lh(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * reshape(F_hat, [1, 3 * nf]);
Lh(:, (-2 : 0)' + 3 * reLIDs) = Lh(:, (-2 : 0)' + 3 * reLIDs) - bs.phitw_face{1, 1} * reshape(F_hat, [1, 3 * nf]);

%**************************************************************************   
%                source term contributions for DG residual
%**************************************************************************
src = getSourceTerm(msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points, t + bt.c(stage.s) * dt, tc);
Lh = Lh - bs.phitw{1} * reshape(msh.elemJac(:, msh.LElems) .* src, [bs.nep, 3 * msh.nLElems]);

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
