% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, t, dt, stage, M, A, IME, tc, quad1, bs, bt, inviscidFlux, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);
W = A * U(1 : end / 3);
R = A * W;
U = reshape(U, [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));
W = reshape(W, [bs.nb, msh.nLElems]);
R = reshape(R, [bs.nb, msh.nLElems]);

% Evaluate physical quantities 
u = bs.phi{1} * U; 
v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);
w = bs.phi{1} * W;
T = getTemperature(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, tc);

faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
wl = bs.phi_face{1, 2} * W(:, leLIDs);
wr = bs.phi_face{1, 1} * W(:, reLIDs);
vl = ul(2, :) ./ ul(1, :);
vr = ur(2, :) ./ ur(1, :);
Tl = getTemperature(ul(1, :), vl, ul(3, :), wl, tc);
Tr = getTemperature(ur(1, :), vr, ur(3, :), wr, tc);

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

% Evaluate viscous stress tensor
z    = bs.phi{1} * Z;
zl   = bs.phi_face{1, 2} * Z(:, leLIDs); 
zr   = bs.phi_face{1, 1} * Z(:, reLIDs); 
tau  = 4 / 3 / tc.Re * z;
taul = 4 / 3 / tc.Re * zl;
taur = 4 / 3 / tc.Re * zr;

%**************************************************************************
%                           Compute heat flux
%************************************************************************** 
% element contributions
Q = -bs.phitw{2} * T;

% internal face contributions
T_hat = 0.5 * (Tl + Tr);
Q(:, leLIDs) = Q(:, leLIDs) + bs.phitw_face{1, 2} * T_hat;
Q(:, reLIDs) = Q(:, reLIDs) - bs.phitw_face{1, 1} * T_hat;

% Take care of mass mtrix
Q = (-8 * tc.Cv / (3 * tc.We * tc.Pr) * IME * Q) ./ msh.elemJac(:, msh.LElems);

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
    r  = bs.phi{1} * R;
    p  = getPressure(u(:, 1 : 3 : end), T);
    xi = getKortewegStressTensor(u(:, 1 : 3 : end), w, r, tc);
    jE = getEnergyFlux(u(:, 1 : 3 : end), w, z, tc);
    F  = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    F(bs.nep + 1 : 2 * bs.nep, :) = F(bs.nep + 1 : 2 * bs.nep, :) - tau - xi;
    F(2 * bs.nep + 1 : end, :)    = F(2 * bs.nep + 1 : end, :)    - (tau + xi) .* v + q + jE;
    
    Lh = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************
rl  = bs.phi_face{1, 2} * R(:, leLIDs); 
rr  = bs.phi_face{1, 1} * R(:, reLIDs);
pl  = getPressure(ul(1, :), Tl);
pr  = getPressure(ur(1, :), Tr);
xil = getKortewegStressTensor(ul(1, :), wl, rl, tc);
xir = getKortewegStressTensor(ur(1, :), wr, rr, tc);
jEl = getEnergyFlux(ul(1, :), wl, zl, tc);
jEr = getEnergyFlux(ur(1, :), wr, zr, tc);
Fl  = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr  = getF(ur(1, :), ur(2, :), ur(3, :), pr);
Fl(2, :) = Fl(2, :) - taul - xil;
Fr(2, :) = Fr(2, :) - taur - xir;
Fl(3, :) = Fl(3, :) - (taul + xil) .* vl + ql + jEl;
Fr(3, :) = Fr(3, :) - (taur + xir) .* vr + qr + jEr;

% dissipation coefficient 
alpha = getDissCoe(ul(1, :), ur(1, :), vl, vr, Tl, Tr, pl, pr, tc, inviscidFlux);
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
LStage(:, stage.Ls) = reshape(Lh(:, (1 : 3 : end)' + (0 : 2)), [3 * msh.nLElems * bs.nb, 1]);

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
    dp = getPressureDeriv(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, tc);
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 3
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - cdt * reshape(bs.phi{1}' * (dp(:, :, i) .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end
    tmp = bs.phi{1}' * (dp(:, :, 4) .* Lamb(:, msh.nLElems + 1 : end));
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * (A' * tmp(:));

    g = KKT.Umin - [u(:, 1 : 3 : end), p];
    res(3 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end

%**************************************************************************
%                           subroutine   
%************************************************************************** 
% Compute dissipation coefficient
function alpha = getDissCoe(rhol, rhor, vl, vr, Tl, Tr, pl, pr, tc, inviscidFlux)

alphal = getEigenMax(rhol, vl, Tl, pl, tc);
alphar = getEigenMax(rhor, vr, Tr, pr, tc);

switch inviscidFlux
    case 1
        alpha = max(max(alphal, alphar));
    case 2
        alpha = max(alphal, alphar);
end

end
