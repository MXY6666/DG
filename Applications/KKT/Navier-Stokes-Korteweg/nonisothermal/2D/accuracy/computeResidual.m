% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, t, dt, stage, M, Ax, Ay, IME, tc, quad2, bs, bt, inviscidFlux, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U  = getStageSolution(U0, KStage, dt, stage, bt);
W1 = Ax * U(1 : end / 4);
W2 = Ay * U(1 : end / 4);
R  = Ax * W1 + Ay * W2;
U  = reshape(U, [bs.nb, 4 * msh.nLElems]);
U  = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));
W1 = reshape(W1, [bs.nb, msh.nLElems]);
W2 = reshape(W2, [bs.nb, msh.nLElems]);
R  = reshape(R, [bs.nb, msh.nLElems]);

% Evaluate physical quantities 
u  = bs.phi{1} * U; 
v  = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];
w1 = bs.phi{1} * W1;
w2 = bs.phi{1} * W2;
T  = getTemperature(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, tc);

% vertical faces
faceIDs = md.intLFaces{1, 3};
leLIDs1 = msh.faceElems(1, faceIDs);
reLIDs1 = msh.faceElems(2, faceIDs);
J1      = msh.faceJac(:, faceIDs);
nf1     = length(faceIDs);

ul1  = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs1), [4 * bs.nfp, nf1]);
ur1  = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs1), [4 * bs.nfp, nf1]);
vl1  = ul1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul1(1 : bs.nfp, :), [2, 1]);
vr1  = ur1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur1(1 : bs.nfp, :), [2, 1]);
w1l1 = bs.phi_face{1, 2} * W1(:, leLIDs1);
w1r1 = bs.phi_face{1, 1} * W1(:, reLIDs1);
w2l1 = bs.phi_face{1, 2} * W2(:, leLIDs1);
w2r1 = bs.phi_face{1, 1} * W2(:, reLIDs1);
Tl1  = getTemperature(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), w1l1, w2l1, tc);
Tr1  = getTemperature(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), w1r1, w2r1, tc);

% horizontal faces
faceIDs = md.intLFaces{2, 3};
leLIDs2 = msh.faceElems(1, faceIDs);
reLIDs2 = msh.faceElems(2, faceIDs);
J2      = msh.faceJac(:, faceIDs);
nf2     = length(faceIDs);

ul2 = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs2), [4 * bs.nfp, nf2]);
ur2 = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs2), [4 * bs.nfp, nf2]);
vl2 = ul2(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul2(1 : bs.nfp, :), [2, 1]);
vr2 = ur2(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur2(1 : bs.nfp, :), [2, 1]);
w1l2 = bs.phi_face{1, 4} * W1(:, leLIDs2);
w1r2 = bs.phi_face{1, 3} * W1(:, reLIDs2);
w2l2 = bs.phi_face{1, 4} * W2(:, leLIDs2);
w2r2 = bs.phi_face{1, 3} * W2(:, reLIDs2);
Tl2  = getTemperature(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), w1l2, w2l2, tc);
Tr2  = getTemperature(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), w1r2, w2r2, tc);

%**************************************************************************
%                     Compute viscous stress tensor
%************************************************************************** 
% element contributions
Z11 = -msh.elemJxix( :, msh.LElems) .* (bs.phitw{2} * v(:, 1 : end / 2));
Z12 = -msh.elemJetay(:, msh.LElems) .* (bs.phitw{3} * v(:, 1 : end / 2));
Z21 = -msh.elemJxix( :, msh.LElems) .* (bs.phitw{2} * v(:, end / 2 + 1 : end));
Z22 = -msh.elemJetay(:, msh.LElems) .* (bs.phitw{3} * v(:, end / 2 + 1 : end));

% internal vertical face contributions
v_hat = 0.5 * (vl1 + vr1);
Z11(:, leLIDs1) = Z11(:, leLIDs1) + J1 .* (bs.phitw_face{1, 2} * v_hat(1 : bs.nfp, :));
Z11(:, reLIDs1) = Z11(:, reLIDs1) - J1 .* (bs.phitw_face{1, 1} * v_hat(1 : bs.nfp, :));
Z21(:, leLIDs1) = Z21(:, leLIDs1) + J1 .* (bs.phitw_face{1, 2} * v_hat(bs.nfp + 1 : end, :));
Z21(:, reLIDs1) = Z21(:, reLIDs1) - J1 .* (bs.phitw_face{1, 1} * v_hat(bs.nfp + 1 : end, :));

% internal horizontal face contributions
v_hat = 0.5 * (vl2 + vr2);
Z12(:, leLIDs2) = Z12(:, leLIDs2) + J2 .* (bs.phitw_face{1, 4} * v_hat(1 : bs.nfp, :));
Z12(:, reLIDs2) = Z12(:, reLIDs2) - J2 .* (bs.phitw_face{1, 3} * v_hat(1 : bs.nfp, :));
Z22(:, leLIDs2) = Z22(:, leLIDs2) + J2 .* (bs.phitw_face{1, 4} * v_hat(bs.nfp + 1 : end, :));
Z22(:, reLIDs2) = Z22(:, reLIDs2) - J2 .* (bs.phitw_face{1, 3} * v_hat(bs.nfp + 1 : end, :));

% Take care of mass mtrix
Z11 = (IME * Z11) ./ msh.elemJac(:, msh.LElems);
Z12 = (IME * Z12) ./ msh.elemJac(:, msh.LElems);
Z21 = (IME * Z21) ./ msh.elemJac(:, msh.LElems);
Z22 = (IME * Z22) ./ msh.elemJac(:, msh.LElems);
TRZ = Z11 + Z22;

% viscous stress tensor
TAU11 = (4 / 3 * Z11 - 2 / 3 * Z22) / tc.Re;
TAU12 = (Z12 +  Z21) / tc.Re;
TAU22 = (-2 / 3 * Z11 + 4 / 3 * Z22) / tc.Re;

%**************************************************************************
%                           Compute heat flux
%************************************************************************** 
% element contributions
Q1 = -msh.elemJxix( :, msh.LElems) .* (bs.phitw{2} * T);
Q2 = -msh.elemJetay(:, msh.LElems) .* (bs.phitw{3} * T);

% internal vertical face contributions
T_hat = 0.5 * (Tl1 + Tr1);
Q1(:, leLIDs1) = Q1(:, leLIDs1) + J1 .* (bs.phitw_face{1, 2} * T_hat);
Q1(:, reLIDs1) = Q1(:, reLIDs1) - J1 .* (bs.phitw_face{1, 1} * T_hat);

% internal horizontal face contributions
T_hat = 0.5 * (Tl2 + Tr2);
Q2(:, leLIDs2) = Q2(:, leLIDs2) + J2 .* (bs.phitw_face{1, 4} * T_hat);
Q2(:, reLIDs2) = Q2(:, reLIDs2) - J2 .* (bs.phitw_face{1, 3} * T_hat);

% Take care of mass mtrix
Q1 = (-8 * tc.Cv / (3 * tc.We * tc.Pr) * IME * Q1) ./ msh.elemJac(:, msh.LElems);
Q2 = (-8 * tc.Cv / (3 * tc.We * tc.Pr) * IME * Q2) ./ msh.elemJac(:, msh.LElems);

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
if (bs.type == 100)
    Lh = zeros(1, 4 * msh.nLElems);
else
    r     = bs.phi{1} * R;
    trz   = bs.phi{1} * TRZ;
    tau11 = bs.phi{1} * TAU11;
    tau12 = bs.phi{1} * TAU12;
    tau22 = bs.phi{1} * TAU22;
    q1    = bs.phi{1} * Q1;
    q2    = bs.phi{1} * Q2;
    p     = getPressure(u(:, 1 : 4 : end), T);
    [xi11, xi12, xi22] = getKortewegStressTensor(u(:, 1 : 4 : end), w1, w2, r, 12, tc);
    [jE1, jE2] = getEnergyFlux(u(:, 1 : 4 : end), w1, w2, trz, 12, tc);
    
    F = getF(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);
    G = getG(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);
    F(bs.nep + 1 : 2 * bs.nep, :)     = F(bs.nep + 1 : 2 * bs.nep, :)     - tau11 - xi11;
    F(2 * bs.nep + 1 : 3 * bs.nep, :) = F(2 * bs.nep + 1 : 3 * bs.nep, :) - tau12 - xi12;
    F(3 * bs.nep + 1 : end, :)        = F(3 * bs.nep + 1 : end, :)        - (tau11 + xi11) .* v(:, 1 : end / 2) - (tau12 + xi12) .* v(:, end / 2 + 1 : end) + q1 + jE1;
    G(bs.nep + 1 : 2 * bs.nep, :)     = G(bs.nep + 1 : 2 * bs.nep, :)     - tau12 - xi12;
    G(2 * bs.nep + 1 : 3 * bs.nep, :) = G(2 * bs.nep + 1 : 3 * bs.nep, :) - tau22 - xi22;
    G(3 * bs.nep + 1 : end, :)        = G(3 * bs.nep + 1 : end, :)        - (tau12 + xi12) .* v(:, 1 : end / 2) - (tau22 + xi22) .* v(:, end / 2 + 1 : end) + q2 + jE2;
    
    Lh = -bs.phitw{2} * reshape(msh.elemJxix( :, msh.LElems) .* F, [bs.nep, 4 * msh.nLElems]) ...                  
         -bs.phitw{3} * reshape(msh.elemJetay(:, msh.LElems) .* G, [bs.nep, 4 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
% internal vertical face contributions
rl     = bs.phi_face{1, 2} * R(:, leLIDs1);
rr     = bs.phi_face{1, 1} * R(:, reLIDs1);
trzl   = bs.phi_face{1, 2} * TRZ(:, leLIDs1);
trzr   = bs.phi_face{1, 1} * TRZ(:, reLIDs1);
tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs1);
tau11r = bs.phi_face{1, 1} * TAU11(:, reLIDs1);
tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs1);
tau12r = bs.phi_face{1, 1} * TAU12(:, reLIDs1);
ql     = bs.phi_face{1, 2} * Q1(:, leLIDs1);
qr     = bs.phi_face{1, 1} * Q1(:, reLIDs1);
pl     = getPressure(ul1(1 : bs.nfp, :), Tl1);
pr     = getPressure(ur1(1 : bs.nfp, :), Tr1);
[xi11l, xi12l] = getKortewegStressTensor(ul1(1 : bs.nfp, :), w1l1, w2l1, rl, 1, tc);
[xi11r, xi12r] = getKortewegStressTensor(ur1(1 : bs.nfp, :), w1r1, w2r1, rr, 1, tc);
jEl = getEnergyFlux(ul1(1 : bs.nfp, :), w1l1, [], trzl, 1, tc);
jEr = getEnergyFlux(ur1(1 : bs.nfp, :), w1r1, [], trzr, 1, tc);

Fl = getF(ul1(1 : bs.nfp, :), ul1(bs.nfp + 1 : 2 * bs.nfp, :), ul1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul1(3 * bs.nfp + 1 : end, :), pl); 
Fr = getF(ur1(1 : bs.nfp, :), ur1(bs.nfp + 1 : 2 * bs.nfp, :), ur1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur1(3 * bs.nfp + 1 : end, :), pr);
Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11l - xi11l;
Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11r - xi11r;
Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12l - xi12l;
Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12r - xi12r;
Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - (tau11l + xi11l) .* vl1(1 : bs.nfp, :) - (tau12l + xi12l) .* vl1(bs.nfp + 1 : end, :) + ql + jEl;
Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - (tau11r + xi11r) .* vr1(1 : bs.nfp, :) - (tau12r + xi12r) .* vr1(bs.nfp + 1 : end, :) + qr + jEr;

% dissipation coefficient 
alpha = getDissCoe(ul1(1 : bs.nfp, :), ur1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), Tl1, Tr1, pl, pr, tc, inviscidFlux);
waveVel = max(alpha, [], 'all');

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur1 - ul1));
Lh(:, (-3 : 0)' + 4 * leLIDs1) = Lh(:, (-3 : 0)' + 4 * leLIDs1) + bs.phitw_face{1, 2} * reshape(J1 .* F_hat, [bs.nfp, 4 * nf1]);
Lh(:, (-3 : 0)' + 4 * reLIDs1) = Lh(:, (-3 : 0)' + 4 * reLIDs1) - bs.phitw_face{1, 1} * reshape(J1 .* F_hat, [bs.nfp, 4 * nf1]);

% internal horizontal face contributions
rl     = bs.phi_face{1, 4} * R(:, leLIDs2);
rr     = bs.phi_face{1, 3} * R(:, reLIDs2);
trzl   = bs.phi_face{1, 4} * TRZ(:, leLIDs2);
trzr   = bs.phi_face{1, 3} * TRZ(:, reLIDs2);
tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs2);
tau12r = bs.phi_face{1, 3} * TAU12(:, reLIDs2);
tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs2);
tau22r = bs.phi_face{1, 3} * TAU22(:, reLIDs2);
ql     = bs.phi_face{1, 4} * Q2(:, leLIDs2);
qr     = bs.phi_face{1, 3} * Q2(:, reLIDs2);
pl     = getPressure(ul2(1 : bs.nfp, :), Tl2);
pr     = getPressure(ur2(1 : bs.nfp, :), Tr2);
[~, xi12l, xi22l] = getKortewegStressTensor(ul2(1 : bs.nfp, :), w1l2, w2l2, rl, 2, tc);
[~, xi12r, xi22r] = getKortewegStressTensor(ur2(1 : bs.nfp, :), w1r2, w2r2, rr, 2, tc);
[~, jEl] = getEnergyFlux(ul2(1 : bs.nfp, :), [], w2l2, trzl, 2, tc);
[~, jEr] = getEnergyFlux(ur2(1 : bs.nfp, :), [], w2r2, trzr, 2, tc);

Fl = getG(ul2(1 : bs.nfp, :), ul2(bs.nfp + 1 : 2 * bs.nfp, :), ul2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul2(3 * bs.nfp + 1 : end, :), pl); 
Fr = getG(ur2(1 : bs.nfp, :), ur2(bs.nfp + 1 : 2 * bs.nfp, :), ur2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur2(3 * bs.nfp + 1 : end, :), pr);
Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12l - xi12l;
Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12r - xi12r;
Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22l - xi22l;
Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22r - xi22r;
Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - (tau12l + xi12l) .* vl2(1 : bs.nfp, :) - (tau22l + xi22l) .* vl2(bs.nfp + 1 : end, :) + ql + jEl;
Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - (tau12r + xi12r) .* vr2(1 : bs.nfp, :) - (tau22r + xi22r) .* vr2(bs.nfp + 1 : end, :) + qr + jEr;

% dissipation coefficient 
alpha = getDissCoe(ul2(1 : bs.nfp, :), ur2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), vr2(bs.nfp + 1 : end, :), Tl2, Tr2, pl, pr, tc, inviscidFlux);
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur2 - ul2));
Lh(:, (-3 : 0)' + 4 * leLIDs2) = Lh(:, (-3 : 0)' + 4 * leLIDs2) + bs.phitw_face{1, 4} * reshape(J2 .* F_hat, [bs.nfp, 4 * nf2]);
Lh(:, (-3 : 0)' + 4 * reLIDs2) = Lh(:, (-3 : 0)' + 4 * reLIDs2) - bs.phitw_face{1, 3} * reshape(J2 .* F_hat, [bs.nfp, 4 * nf2]);

%**************************************************************************   
%                source term contributions for DG residual
%**************************************************************************
ct  = msh.elemCenter(:, msh.LElems);
h   = msh.elemLength(:, msh.LElems);
src = getSourceTerm(ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2), t + bt.c(stage.s) * dt, tc);
Lh = Lh - bs.phitw{1} * reshape(msh.elemJac(:, msh.LElems) .* src, [bs.nep, 4 * msh.nLElems]);

%**************************************************************************   
%                          residual vector
%************************************************************************** 
% Compute DG residual for current Runge-Kutta stage
LStage(:, stage.Ls) = reshape(Lh(:, (1 : 4 : end)' + (0 : 3)), [4 * msh.nLElems * bs.nb, 1]);

% DG discretization contribution for residual vector
res = zeros(4 * KKT.NDOF + 2 * KKT.LDOF, 1);
if (bt.form == 1) 
    res(1 : 4 * KKT.NDOF) = M * (KStage(:, 1) - U0) + dt * (LStage(:, 1 : stage.s) * bt.A(stage.s, 1 : stage.s)');
elseif (bt.form == 2)  
    res(1 : 4 * KKT.NDOF) = M * KStage(:, stage.s) + LStage(:, 1);
end

% Add positivity constraints contribution
if (KKT.positivity)
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    if (bt.form == 1) 
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;       
    end
    dp = getPressureDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, tc);
    
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 4
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - cdt * reshape(bs.phi{1}' * (dp(:, :, i) .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end
    tmp = bs.phi{1}' * (dp(:, :, 5) .* Lamb(:, msh.nLElems + 1 : end));
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * (Ax' * tmp(:));
    tmp = bs.phi{1}' * (dp(:, :, 6) .* Lamb(:, msh.nLElems + 1 : end));
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * (Ay' * tmp(:));

    g = KKT.Umin - [u(:, 1 : 4 : end), p];
    res(4 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end

%**************************************************************************
%                           subroutine   
%************************************************************************** 
% Compute dissipation coefficient
function alpha = getDissCoe(rhol, rhor, vnl, vnr, Tl, Tr, pl, pr, tc, inviscidFlux)

alphal = getEigenMax(rhol, vnl, Tl, pl, tc);
alphar = getEigenMax(rhor, vnr, Tr, pr, tc);

switch inviscidFlux
    case 1
        alpha = repmat(max(max(alphal, alphar), [], 'all'), [size(rhol, 1), 1]);
    case 2
        alpha = max(alphal, alphar);
end

end
