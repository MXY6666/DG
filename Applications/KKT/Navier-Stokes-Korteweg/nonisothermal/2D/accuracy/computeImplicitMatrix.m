function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, Ax, Ay, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, bs, bt, inviscidFlux, KKT, newton, mtol)

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
u = bs.phi{1} * U; 
v = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];

% vertical faces
faceIDs = md.intLFaces{1, 3};
leIDs1  = msh.faceElems(1, faceIDs);
reIDs1  = msh.faceElems(2, faceIDs);
leLIDs1 = leIDs1;
reLIDs1 = reIDs1;
J1      = msh.faceJac(:, faceIDs);
nf1     = length(faceIDs);

ul1 = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs1), [4 * bs.nfp, nf1]);
ur1 = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs1), [4 * bs.nfp, nf1]);
vl1 = ul1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul1(1 : bs.nfp, :), [2, 1]);
vr1 = ur1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur1(1 : bs.nfp, :), [2, 1]);

% horizontal faces
faceIDs = md.intLFaces{2, 3};
leIDs2  = msh.faceElems(1, faceIDs);
reIDs2  = msh.faceElems(2, faceIDs);
leLIDs2 = leIDs2;
reLIDs2 = reIDs2;
J2      = msh.faceJac(:, faceIDs);
nf2     = length(faceIDs);

ul2 = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs2), [4 * bs.nfp, nf2]);
ur2 = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs2), [4 * bs.nfp, nf2]);
vl2 = ul2(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul2(1 : bs.nfp, :), [2, 1]);
vr2 = ur2(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur2(1 : bs.nfp, :), [2, 1]);

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);

r11 = zeros((16 * (msh.nLElems + 4 * md.nIntLFaces) + 4 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

r12 = zeros(length(pr12) + 15 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c12 = r12;
v12 = r12;

r21 = zeros(length(pr21) + 16 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c21 = r21;
v21 = r21;

r22i = zeros(length(pr22i) + 4 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c22i = r22i;
v22i = r22i;

e12 = length(pr12);
r12(1 : e12) = pr12;
c12(1 : e12) = pc12;
v12(1 : e12) = pv12;

e21 = length(pr21);
r21(1 : e21) = pr21;
c21(1 : e21) = pc21;
v21(1 : e21) = pv21;

e22 = length(pr22i);
r22i(1 : e22) = pr22i;
c22i(1 : e22) = pc22i;
v22i(1 : e22) = pv22i;

%**************************************************************************
%                  Compute viscous stress tensor
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
%                       element contributions   
%**************************************************************************
w1    = bs.phi{1} * W1;
w2    = bs.phi{1} * W2;
r     = bs.phi{1} * R;
trz   = bs.phi{1} * TRZ;
tau11 = bs.phi{1} * TAU11;
tau12 = bs.phi{1} * TAU12;
tau22 = bs.phi{1} * TAU22;
T     = getTemperature(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, tc);
p     = getPressure(u(:, 1 : 4 : end), T);
[xi11, xi12, xi22] = getKortewegStressTensor(u(:, 1 : 4 : end), w1, w2, r, 12, tc);
dT = getTemperatureDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, tc);
dp = getPressureDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, tc);

% block-(1,1)
e11 = 0;
for n = 1 : 16
    if (n == 1 || n == 4)
        continue;
    end
    i = ceil(n / 4); j = n - (i - 1) * 4;
    ent1 = getJacobianEntry1(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), p, w1, w2, r, trz, tau11, tau12, xi11, xi12, dp, i, j, tc);
    ent2 = getJacobianEntry2(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), p, w1, w2, r, trz, tau12, tau22, xi12, xi22, dp, i, j, tc);
    
    s11 = e11 + 1; l = msh.nLElems * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s11 : e11) = -bs.phitw{2} * (repelem(msh.elemJxix( :, msh.LElems) .* ent1, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])) + ...
                     -bs.phitw{3} * (repelem(msh.elemJetay(:, msh.LElems) .* ent2, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
end

% block-(1,2)
for n = 1 : 36
    i = ceil(n / 9); j = n - (i - 1) * 9;
    [ent1, flag1] = getJacobianEntry1(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), p, w1, w2, r, trz, tau11, tau12, xi11, xi12, dp, i, j + 4, tc);
    [ent2, flag2] = getJacobianEntry2(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), p, w1, w2, r, trz, tau12, tau22, xi12, xi22, dp, i, j + 4, tc);
    
    if (flag1 == 2 || flag2 == 2)
        s12 = e12 + 1; l = msh.nLElems * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        v12(s12 : e12) = -bs.phitw{2} * (repelem(msh.elemJxix( :, msh.LElems) .* ent1, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])) + ...
                         -bs.phitw{3} * (repelem(msh.elemJetay(:, msh.LElems) .* ent2, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    end
end

% block-(2,1)
for n = 1 : 12
    i = ceil(n / 4); j = n - (i - 1) * 4;
    i1 = 2 * i + 2; i2 = i1 + 1; 
    if (i <= 2)
        [ent, flag] = getVelocityDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), i, j);
    elseif (i == 3)
        ent = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dT(:, :, j);
        flag = 2;
    end
    if (flag == 0)
       continue; 
    end    
    
    s21 = e21 + 1; l = msh.nLElems * bs.nb2; e21 = e21 + l;    
    r21(s21 : e21) = re1 + ((i1 - 1) * msh.nLElems : i1 * msh.nLElems - 1) * bs.nb;
    c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems : j  * msh.nLElems - 1) * bs.nb;
    v21(s21 : e21) = reshape(-2 * IME * bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemLength(1, msh.LElems); 
    
    s21 = e21 + 1; l = msh.nLElems * bs.nb2; e21 = e21 + l; 
    r21(s21 : e21) = re1 + ((i2 - 1) * msh.nLElems : i2 * msh.nLElems - 1) * bs.nb;
    c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems : j  * msh.nLElems - 1) * bs.nb;
    v21(s21 : e21) = reshape(-2 * IME * bs.phitw{3} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemLength(2, msh.LElems);     
end

% inverse of block-(2,2)
for j = 1 : 2
    i1 = 8; i2 = 9;
    ent = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dT(:, :, j + 4);
    
    s22 = e22 + 1; l = msh.nLElems * bs.nb2; e22 = e22 + l;  
    r22i(s22 : e22) = re1 + ((i1 - 1) * msh.nLElems : i1 * msh.nLElems - 1) * bs.nb;
    c22i(s22 : e22) = ce1 + ((j  - 1) * msh.nLElems : j  * msh.nLElems - 1) * bs.nb;
    v22i(s22 : e22) = reshape(2 * IME * bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemLength(1, msh.LElems);
    
    s22 = e22 + 1; l = msh.nLElems * bs.nb2; e22 = e22 + l;  
    r22i(s22 : e22) = re1 + ((i2 - 1) * msh.nLElems : i2 * msh.nLElems - 1) * bs.nb;
    c22i(s22 : e22) = ce1 + ((j  - 1) * msh.nLElems : j  * msh.nLElems - 1) * bs.nb;
    v22i(s22 : e22) = reshape(2 * IME * bs.phitw{3} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemLength(2, msh.LElems);    
end

%**************************************************************************
%                        face contributions   
%**************************************************************************
% internal vertical face contributions
w1l    = bs.phi_face{1, 2} * W1(:, leLIDs1);
w1r    = bs.phi_face{1, 1} * W1(:, reLIDs1);
w2l    = bs.phi_face{1, 2} * W2(:, leLIDs1);
w2r    = bs.phi_face{1, 1} * W2(:, reLIDs1);
rl     = bs.phi_face{1, 2} * R(:, leLIDs1);
rr     = bs.phi_face{1, 1} * R(:, reLIDs1);
trzl   = bs.phi_face{1, 2} * TRZ(:, leLIDs1);
trzr   = bs.phi_face{1, 1} * TRZ(:, reLIDs1);
tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs1);
tau11r = bs.phi_face{1, 1} * TAU11(:, reLIDs1);
tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs1);
tau12r = bs.phi_face{1, 1} * TAU12(:, reLIDs1);
[xi11l, xi12l] = getKortewegStressTensor(ul1(1 : bs.nfp, :), w1l, w2l, rl, 1, tc);
[xi11r, xi12r] = getKortewegStressTensor(ur1(1 : bs.nfp, :), w1r, w2r, rr, 1, tc);
Tl  = getTemperature(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
Tr  = getTemperature(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);
pl  = getPressure(ul1(1 : bs.nfp, :), Tl);
pr  = getPressure(ur1(1 : bs.nfp, :), Tr);
dTl = getTemperatureDeriv(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
dTr = getTemperatureDeriv(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);
dpl = getPressureDeriv(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
dpr = getPressureDeriv(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul1(1 : bs.nfp, :), ur1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), Tl, Tr, pl, pr, dTl, dTr, dpl, dpr, 1, tc, inviscidFlux);
    
% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;    
    entl = getJacobianEntry1(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), pl, w1l, w2l, rl, trzl, tau11l, tau12l, xi11l, xi12l, dpl, i, j, tc);
    entr = getJacobianEntry1(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), pr, w1r, w2r, rr, trzr, tau11r, tau12r, xi11r, xi12r, dpr, i, j, tc);   
    H_ul = 0.5 * (entl - dalphal(:, :, j) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar(:, :, j) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s11 = e11 + 1; l = 4 * nf1 * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;  
    v11(s11 : e11) = [bs.phitw_face{1, 2} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                      bs.phitw_face{1, 2} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];    
end

% block-(1,2)
for n = 1 : 36
    i = ceil(n / 9); j = n - (i - 1) * 9;   
    entl         = getJacobianEntry1(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), pl, w1l, w2l, rl, trzl, tau11l, tau12l, xi11l, xi12l, dpl, i, j + 4, tc);
    [entr, flag] = getJacobianEntry1(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), pr, w1r, w2r, rr, trzr, tau11r, tau12r, xi11r, xi12r, dpr, i, j + 4, tc);  
    
    if (j == 1 || j == 2)
        H_ul = 0.5 * (entl - dalphal(:, :, j + 4) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)));
        H_ur = 0.5 * (entr - dalphar(:, :, j + 4) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)));
        
        s12 = e12 + 1; l = 4 * nf1 * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
        c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;  
        v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                          bs.phitw_face{1, 2} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                         -bs.phitw_face{1, 1} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                         -bs.phitw_face{1, 1} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];         
    else
        if (flag == 2)
            H_ul = 0.5 * entl;
            H_ur = 0.5 * entr;

            s12 = e12 + 1; l = 4 * nf1 * bs.nb2; e12 = e12 + l;
            r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
            c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;  
            v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                              bs.phitw_face{1, 2} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                             -bs.phitw_face{1, 1} * (repelem(J1 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                             -bs.phitw_face{1, 1} * (repelem(J1 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];  
        end
    end
end

% block-(2,1)
for n = 1 : 12
    i = ceil(n / 4); j = n - (i - 1) * 4;
    i1 = 2 * i + 2; 
    if (i <= 2)
        entl         = getVelocityDeriv(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), i, j);
        [entr, flag] = getVelocityDeriv(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), i, j);
    elseif (i == 3)
        entl = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, j);
        entr = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, j);
        flag = 2;
    end
    if (flag == 0)
       continue; 
    end  
    
    s21 = e21 + 1; l = 4 * nf1 * bs.nb2; e21 = e21 + l;    
    r21(s21 : e21) = re1 + [(i1 - 1) * msh.nLElems + leLIDs1 - 1, (i1 - 1) * msh.nLElems + leLIDs1 - 1, (i1 - 1) * msh.nLElems + reLIDs1 - 1, (i1 - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c21(s21 : e21) = ce1 + [(j  - 1) * msh.nLElems + leLIDs1 - 1, (j  - 1) * msh.nLElems + reLIDs1 - 1, (j  - 1) * msh.nLElems + leLIDs1 - 1, (j  - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    v21(s21 : e21) = [reshape( IME * bs.phitw_face{1, 2} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, leIDs1), ...
                      reshape( IME * bs.phitw_face{1, 2} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, leIDs1), ...
                      reshape(-IME * bs.phitw_face{1, 1} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, reIDs1), ...
                      reshape(-IME * bs.phitw_face{1, 1} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, reIDs1)];
end

% inverse of block-(2,2)
i = 8;
for j = 1 : 2
    entl = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, j + 4);
    entr = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, j + 4);
    
    s22 = e22 + 1; l = 4 * nf1 * bs.nb2; e22 = e22 + l;    
    r22i(s22 : e22) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c22i(s22 : e22) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    v22i(s22 : e22) = [reshape(-IME * bs.phitw_face{1, 2} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, leIDs1), ...
                       reshape(-IME * bs.phitw_face{1, 2} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, leIDs1), ...
                       reshape( IME * bs.phitw_face{1, 1} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, reIDs1), ...
                       reshape( IME * bs.phitw_face{1, 1} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(1, reIDs1)];    
end

% internal horizontal face contributions
w1l    = bs.phi_face{1, 4} * W1(:, leLIDs2);
w1r    = bs.phi_face{1, 3} * W1(:, reLIDs2);
w2l    = bs.phi_face{1, 4} * W2(:, leLIDs2);
w2r    = bs.phi_face{1, 3} * W2(:, reLIDs2);
rl     = bs.phi_face{1, 4} * R(:, leLIDs2);
rr     = bs.phi_face{1, 3} * R(:, reLIDs2);
trzl   = bs.phi_face{1, 4} * TRZ(:, leLIDs2);
trzr   = bs.phi_face{1, 3} * TRZ(:, reLIDs2);
tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs2);
tau12r = bs.phi_face{1, 3} * TAU12(:, reLIDs2);
tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs2);
tau22r = bs.phi_face{1, 3} * TAU22(:, reLIDs2);
[~, xi12l, xi22l] = getKortewegStressTensor(ul2(1 : bs.nfp, :), w1l, w2l, rl, 2, tc);
[~, xi12r, xi22r] = getKortewegStressTensor(ur2(1 : bs.nfp, :), w1r, w2r, rr, 2, tc);
Tl  = getTemperature(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
Tr  = getTemperature(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);
pl  = getPressure(ul2(1 : bs.nfp, :), Tl);
pr  = getPressure(ur2(1 : bs.nfp, :), Tr);
dTl = getTemperatureDeriv(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
dTr = getTemperatureDeriv(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);
dpl = getPressureDeriv(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), w1l, w2l, tc);
dpr = getPressureDeriv(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), w1r, w2r, tc);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul2(1 : bs.nfp, :), ur2(1 : bs.nfp, :), vl2(bs.nfp + 1 : 2 * bs.nfp, :), vr2(bs.nfp + 1 : 2 * bs.nfp, :), Tl, Tr, pl, pr, dTl, dTr, dpl, dpr, 2, tc, inviscidFlux);
              
% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4; 
    entl = getJacobianEntry2(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), pl, w1l, w2l, rl, trzl, tau12l, tau22l, xi12l, xi22l, dpl, i, j, tc);
    entr = getJacobianEntry2(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), pr, w1r, w2r, rr, trzr, tau12r, tau22r, xi12r, xi22r, dpr, i, j, tc);   
    H_ul = 0.5 * (entl - dalphal(:, :, j) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar(:, :, j) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s11 = e11 + 1; l = 4 * nf2 * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;  
    v11(s11 : e11) = [bs.phitw_face{1, 4} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                      bs.phitw_face{1, 4} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), ...
                     -bs.phitw_face{1, 3} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                     -bs.phitw_face{1, 3} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2]))];    
end

% block-(1,2)
for n = 1 : 36
    i = ceil(n / 9); j = n - (i - 1) * 9; 
    entl         = getJacobianEntry2(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), pl, w1l, w2l, rl, trzl, tau12l, tau22l, xi12l, xi22l, dpl, i, j + 4, tc);
    [entr, flag] = getJacobianEntry2(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), pr, w1r, w2r, rr, trzr, tau12r, tau22r, xi12r, xi22r, dpr, i, j + 4, tc);  
    
    if (j == 1 || j == 2)
        H_ul = 0.5 * (entl - dalphal(:, :, j + 4) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)));
        H_ur = 0.5 * (entr - dalphar(:, :, j + 4) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)));
        
        s12 = e12 + 1; l = 4 * nf2 * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
        c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;  
        v12(s12 : e12) = [bs.phitw_face{1, 4} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                          bs.phitw_face{1, 4} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), ...
                         -bs.phitw_face{1, 3} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                         -bs.phitw_face{1, 3} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2]))];         
    else
        if (flag == 2)
            H_ul = 0.5 * entl;
            H_ur = 0.5 * entr;

            s12 = e12 + 1; l = 4 * nf2 * bs.nb2; e12 = e12 + l;
            r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
            c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;  
            v12(s12 : e12) = [bs.phitw_face{1, 4} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                              bs.phitw_face{1, 4} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), ...
                             -bs.phitw_face{1, 3} * (repelem(J2 .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                             -bs.phitw_face{1, 3} * (repelem(J2 .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2]))];   
        end
    end
end

% block-(2,1)
for n = 1 : 12
    i = ceil(n / 4); j = n - (i - 1) * 4;
    i2 = 2 * i + 3; 
    if (i <= 2)
        entl         = getVelocityDeriv(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), i, j);
        [entr, flag] = getVelocityDeriv(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), i, j);
    elseif (i == 3)
        entl = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, j);
        entr = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, j);
        flag = 2;
    end
    if (flag == 0)
       continue; 
    end  
    
    s21 = e21 + 1; l = 4 * nf2 * bs.nb2; e21 = e21 + l;    
    r21(s21 : e21) = re1 + [(i2 - 1) * msh.nLElems + leLIDs2 - 1, (i2 - 1) * msh.nLElems + leLIDs2 - 1, (i2 - 1) * msh.nLElems + reLIDs2 - 1, (i2 - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    c21(s21 : e21) = ce1 + [(j  - 1) * msh.nLElems + leLIDs2 - 1, (j  - 1) * msh.nLElems + reLIDs2 - 1, (j  - 1) * msh.nLElems + leLIDs2 - 1, (j  - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    v21(s21 : e21) = [reshape( IME * bs.phitw_face{1, 4} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, leIDs2), ...
                      reshape( IME * bs.phitw_face{1, 4} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, leIDs2), ...
                      reshape(-IME * bs.phitw_face{1, 3} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, reIDs2), ...
                      reshape(-IME * bs.phitw_face{1, 3} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, reIDs2)];
end

% inverse of block-(2,2)
i = 9;
for j = 1 : 2
    entl = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, j + 4);
    entr = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, j + 4);
    
    s22 = e22 + 1; l = 4 * nf2 * bs.nb2; e22 = e22 + l;    
    r22i(s22 : e22) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    c22i(s22 : e22) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    v22i(s22 : e22) = [reshape(-IME * bs.phitw_face{1, 4} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, leIDs2), ...
                       reshape(-IME * bs.phitw_face{1, 4} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, leIDs2), ...
                       reshape( IME * bs.phitw_face{1, 3} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, reIDs2), ...
                       reshape( IME * bs.phitw_face{1, 3} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), [bs.nb2, nf2]) ./ msh.elemLength(2, reIDs2)];    
end

v11(1 : e11) = dt * bt.diagCoe * v11(1 : e11); 
v12(1 : e12) = dt * bt.diagCoe * v12(1 : e12); 
%**************************************************************************
%                        mass matrix contributions   
%**************************************************************************
s11 = e11 + 1; l = 4 * msh.nLElems * bs.nb2; e11 = e11 + l;
r11(s11 : e11) = re1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
c11(s11 : e11) = ce1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
v11(s11 : e11) = repmat(msh.elemJac(:, msh.LElems) .* reshape(bs.phitw{1} * bs.phi{1}, [bs.nb2, 1]), [1, 4]);

%**************************************************************************
%                   global matrix for DG discretization
%**************************************************************************
% Assemble
r11(e11  + 1 : end) = [];
c11(e11  + 1 : end) = [];
v11(e11  + 1 : end) = [];
r12(e12  + 1 : end) = [];
c12(e12  + 1 : end) = [];
v12(e12  + 1 : end) = [];
r21(e21  + 1 : end) = [];
c21(e21  + 1 : end) = [];
v21(e21  + 1 : end) = [];
r22i(e22 + 1 : end) = [];
c22i(e22 + 1 : end) = [];
v22i(e22 + 1 : end) = [];
v11  = mychop(v11, mtol);
v12  = mychop(v12, mtol);
v21  = mychop(v21, mtol);
v22i = mychop(v22i, mtol);
M11  = sparse(r11,  c11,  v11,  4 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
M12  = sparse(r12,  c12,  v12,  4 * msh.nLElems * bs.nb, 9 * msh.nLElems * bs.nb);
M21  = sparse(r21,  c21,  v21,  9 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
M22i = sparse(r22i, c22i, v22i, 9 * msh.nLElems * bs.nb, 9 * msh.nLElems * bs.nb);
Imp.Mat = M11 - M12 * M22i * M21;

%**************************************************************************
%                       KKT constraints contribution 
%**************************************************************************
% Add positivity constraints contribution
if (KKT.positivity)      
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    if (bt.form == 1)
        cdt  = 1;
        cdt2 = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;
        cdt2 = cdt * cdt;
    end  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);    
    
    % pre-calculate some matrices
    rtmp = re1 + (0 : msh.nLElems - 1) * bs.nb;
    ctmp = ce1 + (0 : msh.nLElems - 1) * bs.nb;
    
    dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, 1, 5, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp11 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, 1, 6, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp12 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, 5, 5, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp13 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, 6, 6, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp14 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    Mtmp1 = Mtmp11 * Ax + Mtmp12 * Ay + Ax' * (Mtmp11 + Mtmp13 * Ax) + Ay' * (Mtmp12 + Mtmp14 * Ay);
    [rKKT1, cKKT1, vKKT1] = find(Mtmp1);
    
    rtmp  = re2 + (0 : msh.nLElems - 1) * bs.nep;
    ctmp  = ce2 + (0 : msh.nLElems - 1) * bs.nb;  
    
    vtmp = repelem(dp(:, :, 5), 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    Mtmp21 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.LDOF, KKT.NDOF);
    
    vtmp = repelem(dp(:, :, 6), 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    Mtmp22 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.LDOF, KKT.NDOF);
    
    Mtmp2 = Mtmp21 * Ax + Mtmp22 * Ay;
    [rKKT2, cKKT2, vKKT2] = find(Mtmp2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % Initialize row, column and value array for the KKT contributions
    rKKT = zeros(length(vKKT1) + 2 * length(vKKT2) + 9 * msh.nLElems * bs.nb2 + 10 * msh.nLElems * bs.nb * bs.nep + 2 * KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;    
    
    e = length(vKKT1) + 2 * length(vKKT2);
    rKKT(1 : e) = [rKKT1; 4 * KKT.NDOF + KKT.LDOF + rKKT2; cKKT2];
    cKKT(1 : e) = [cKKT1; cKKT2; 4 * KKT.NDOF + KKT.LDOF + rKKT2];
    vKKT(1 : e) = [vKKT1; vKKT2; -vKKT2];
    
    % block (1,1)
    for n = [1, 2, 3, 4, 7, 8, 13, 15, 19]
        i = ceil(n / 6); j = n - (i - 1) * 6;
        dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), w1, w2, i, j, tc);

        s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
        rKKT(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        cKKT(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
    end    
    
    % KKT-block (3,1)    
    s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
    stmp = s;
    rKKT(s : e) = 4 * KKT.NDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
    cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
    vKKT(s : e) = repmat(cdt * bs.phi{1}(:), [msh.nLElems, 1]);
    for j = 1 : 4      
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 4 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp(:, :, j), 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 5 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    [~, set2] = getActiveSets(u(:, 1 : 4 : end), p, Lamb, dU, dLamb, [rKKT(stmp : etmp); 4 * KKT.NDOF + KKT.LDOF + rKKT2], [cKKT(stmp : etmp); cKKT2], [vKKT(stmp : etmp); vKKT2], KKT, newton);    
      
    % KKT-block (3,3)
    s = e + 1; l = 2 * KKT.LDOF; e = e + l;
    rKKT(s : e) = 4 * KKT.NDOF + (1 : 2 * KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    [rDG, cDG, vDG] = find(Imp.Mat);
    vKKT = mychop(vKKT, mtol);
    Imp.Mat = sparse([rDG; rKKT], [cDG; cKKT], [vDG; vKKT], 4 * KKT.NDOF + 2 * KKT.LDOF, 4 * KKT.NDOF + 2 * KKT.LDOF);
    Imp.Mat([false(4 * KKT.NDOF, 1); set2], 1 : 4 * KKT.NDOF) = 0;
     
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix
        Imp.GT  = Imp.Mat';
        Imp.Mat = Imp.GT * Imp.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(Imp.Mat, 1));
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        Imp.GT = speye(size(Imp.Mat, 1));
    end
else
    Imp.GT  = speye(4 * KKT.NDOF);   
end

% Apply row-column scaling
if (newton.scaleMat)
    [Imp.Mat, Imp.PR, Imp.PC] = rowcolScale(Imp.Mat, 1.e-2, 5);
else
    si = size(Imp.Mat, 1);
    Imp.PR = speye(si);
    Imp.PC = speye(si);
end

end

%**************************************************************************
%                           subroutine   
%************************************************************************** 
% Compute dissipation coefficient and its derivatives
function [alpha, dalphal, dalphar] = getDissCoeAndDeriv(rhol, rhor, vnl, vnr, Tl, Tr, pl, pr, dTl, dTr, dpl, dpr, dir, tc, inviscidFlux)

% dissipation coefficient at each face
betal = 2 * (pl ./ rhol + 4 * Tl .* (3 + tc.Cv * (2 * rhol - 3)) ./ (tc.Cv * (rhol - 3).^2));
betar = 2 * (pr ./ rhor + 4 * Tr .* (3 + tc.Cv * (2 * rhor - 3)) ./ (tc.Cv * (rhor - 3).^2));
alphal = max(abs(vnl - sqrt(betal)), abs(vnl + sqrt(betal)));
alphar = max(abs(vnr - sqrt(betar)), abs(vnr + sqrt(betar)));
alphalIsChosen = alphal >= alphar;
alpharIsChosen = ~alphalIsChosen;
alpha = alphar;
alpha(alphalIsChosen) = alphal(alphalIsChosen);

if (inviscidFlux == 1)
    % global dissipation coefficient
    [alpha, ind] = max(alpha, [], 'all', 'linear');
    [ind1, ind2] = ind2sub(size(alphal), ind);
    
    % physical quantities evaluated at the point where global dissipation 
    % coefficient is get
    if (alphal(ind1, ind2) >= alphar(ind1, ind2))
        rho  = rhol(ind1, ind2);
        vn   = vnl(ind1, ind2); 
        T    = Tl(ind1, ind2);
        p    = pl(ind1, ind2);
        dT   = reshape(dTl(ind1, ind2, :), [6, 1]);
        dp   = reshape(dpl(ind1, ind2, :), [6, 1]); 
        beta = betal(ind1, ind2);
    else
        rho  = rhor(ind1, ind2);
        vn   = vnr(ind1, ind2); 
        T    = Tr(ind1, ind2);
        p    = pr(ind1, ind2);
        dT   = reshape(dTr(ind1, ind2, :), [6, 1]);
        dp   = reshape(dpr(ind1, ind2, :), [6, 1]); 
        beta = betar(ind1, ind2);           
    end
    
    % derivative of beta
    tmp = 4 * (3 + tc.Cv * (2 * rho - 3)) ./ (tc.Cv * (rho - 3).^2);
    dbeta = 2 * (dp ./ rho + tmp .* dT);
    dbeta(1) = dbeta(1) - 2 * (p ./ rho.^2 + 8 / tc.Cv * T .* (tc.Cv * rho + 3) ./ (rho - 3).^3);
    
    % derivative of global dissipation coefficient with respect to left and right state
    if (dir == 1)
        I = 2;
    elseif (dir == 2)
        I = 3;
    end
     
    if (beta >= 0)
        sig = sign(vn); sig(sig == 0) = 1;
        dalpha = 0.5 ./ sqrt(beta) .* dbeta;
        dalpha([1, I]) = dalpha([1, I]) + sig ./ rho .* [-vn; 1];
    else
        tmp = 0.5 ./ sqrt(vn.^2 - beta);
        dalpha = -tmp .* dbeta;
        dalpha([1, I]) = dalpha([1, I]) + tmp .* 2 .* vn ./ rho .* [-vn; 1];
    end
    if (alphal(ind1, ind2) >= alphar(ind1, ind2))
        dalphal = reshape(dalpha, [1, 1, 6]);
        dalphar = zeros(1, 1, 6);
    else     
        dalphal = zeros(1, 1, 6);
        dalphar = reshape(dalpha, [1, 1, 6]);
    end  
elseif (inviscidFlux == 2)    
    % derivative of beta
    tmpl = 4 * (3 + tc.Cv * (2 * rhol - 3)) ./ (tc.Cv * (rhol - 3).^2);
    tmpr = 4 * (3 + tc.Cv * (2 * rhor - 3)) ./ (tc.Cv * (rhor - 3).^2);
    dbetal = 2 * (dpl ./ rhol + tmpl .* dTl);
    dbetar = 2 * (dpr ./ rhor + tmpr .* dTr);
    dbetal(:, :, 1) = dbetal(:, :, 1) - 2 * (pl ./ rhol.^2 + 8 / tc.Cv * Tl .* (tc.Cv * rhol + 3) ./ (rhol - 3).^3);
    dbetar(:, :, 1) = dbetar(:, :, 1) - 2 * (pr ./ rhor.^2 + 8 / tc.Cv * Tr .* (tc.Cv * rhor + 3) ./ (rhor - 3).^3);

    % derivative of dissipation coefficient with respect to left and right state
    dalphal = zeros(size(dTl));
    dalphar = zeros(size(dTr));

    betalIsPositive = betal >= 0;
    betarIsPositive = betar >= 0;
    ind1 = repmat(alphalIsChosen & betalIsPositive, [1, 6]);
    ind2 = repmat(alphalIsChosen & (~betalIsPositive), [1, 6]);
    ind3 = repmat(alpharIsChosen & betarIsPositive, [1, 6]);
    ind4 = repmat(alpharIsChosen & (~betarIsPositive), [1, 6]);

    sigl = sign(vnl); 
    sigr = sign(vnr);  
    sigl(sigl == 0) = 1; 
    sigr(sigr == 0) = 1;
    tmpl = 0.5 ./ sqrt(vnl.^2 - betal);
    tmpr = 0.5 ./ sqrt(vnr.^2 - betar);
    if (dir == 1)
        I = 2;
    elseif (dir == 2)
        I = 3;
    end

    dalphal1 = 0.5 ./ sqrt(betal) .* dbetal;
    dalphar1 = 0.5 ./ sqrt(betar) .* dbetar;
    dalphal2 = -tmpl .* dbetal;
    dalphar2 = -tmpr .* dbetar;

    tmp1 = sigl ./ rhol;
    tmp2 = 2 * tmpl .* vnl ./ rhol;
    dalphal1(:, :, 1) = dalphal1(:, :, 1) - tmp1 .* vnl;
    dalphal1(:, :, I) = dalphal1(:, :, I) + tmp1;
    dalphal2(:, :, 1) = dalphal2(:, :, 1) - tmp2 .* vnl;
    dalphal2(:, :, I) = dalphal2(:, :, I) + tmp2;

    tmp1 = sigr ./ rhor;
    tmp2 = 2 * tmpr .* vnr ./ rhor;
    dalphar1(:, :, 1) = dalphar1(:, :, 1) - tmp1 .* vnr;
    dalphar1(:, :, I) = dalphar1(:, :, I) + tmp1;
    dalphar2(:, :, 1) = dalphar2(:, :, 1) - tmp2 .* vnr;
    dalphar2(:, :, I) = dalphar2(:, :, I) + tmp2;

    dalphal(ind1) = dalphal1(ind1);
    dalphal(ind2) = dalphal2(ind2);
    dalphar(ind3) = dalphar1(ind3);
    dalphar(ind4) = dalphar2(ind4);
end

end
