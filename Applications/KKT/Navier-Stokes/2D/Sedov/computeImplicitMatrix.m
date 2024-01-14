function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, isViscous, pr12, pc12, pv12, IME, tc, bs, bt, KKT, newton, mtol)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

% Evaluate physical quantities 
u = bs.phi{1} * U; 
v = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];

% internal vertical faces
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

% internal horizontal faces
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

nBndLFaces = sum(md.nBndLFaces);
r11 = zeros((16 * (msh.nLElems + 4 * md.nIntLFaces + nBndLFaces) + 4 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

if isViscous
    r12 = zeros(length(pr12) + 4 * (msh.nLElems + 4 * md.nIntLFaces + nBndLFaces) * bs.nb2, 1);
    c12 = r12;
    v12 = r12;

    r21 = zeros(16 * (msh.nLElems + 4 * md.nIntLFaces + nBndLFaces) * bs.nb2, 1);
    c21 = r21;
    v21 = r21;

    e12 = length(pr12);
    r12(1 : e12) = pr12;
    c12(1 : e12) = pc12;
    v12(1 : e12) = pv12;
end

%**************************************************************************
%                  Compute viscous stress tensor
%************************************************************************** 
if isViscous
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

    % supersonic outflow boundary face contributions
    faceIDs = md.bndLFaces{2, 2};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul = bs.phi_face{1, 2} * U(:, (-3 : -1)' + 4 * leLIDs);
    vl(:, 2 : 2 : 2 * nf) = ul(:, 3 : 3 : end) ./ ul(:, 1 : 3 : end);
    vl(:, 1 : 2 : end)    = ul(:, 2 : 3 : end) ./ ul(:, 1 : 3 : end);

    Z11(:, leLIDs) = Z11(:, leLIDs) + J .* (bs.phitw_face{1, 2} * vl(:, 1 : 2 : end));
    Z21(:, leLIDs) = Z21(:, leLIDs) + J .* (bs.phitw_face{1, 2} * vl(:, 2 : 2 : end));

    faceIDs = md.bndLFaces{4, 2};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul = bs.phi_face{1, 4} * U(:, (-3 : -1)' + 4 * leLIDs);
    vl(:, 2 : 2 : 2 * nf) = ul(:, 3 : 3 : end) ./ ul(:, 1 : 3 : end);
    vl(:, 1 : 2 : end)    = ul(:, 2 : 3 : end) ./ ul(:, 1 : 3 : end);

    Z12(:, leLIDs) = Z12(:, leLIDs) + J .* (bs.phitw_face{1, 4} * vl(:, 1 : 2 : end));
    Z22(:, leLIDs) = Z22(:, leLIDs) + J .* (bs.phitw_face{1, 4} * vl(:, 2 : 2 : end));

    % slip solid wall boundary face contributions
    faceIDs = md.bndLFaces{1, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);

    ul = bs.phi_face{1, 1} * U(:, [-3; -1] + 4 * leLIDs);
    vl = ul(:, 2 : 2 : end) ./ ul(:, 1 : 2 : end);

    Z21(:, leLIDs) = Z21(:, leLIDs) - J .* (bs.phitw_face{1, 1} * vl);

    faceIDs = md.bndLFaces{3, 1};
    leLIDs  = msh.faceElems(1, faceIDs); 
    J       = msh.faceJac(:, faceIDs);

    ul = bs.phi_face{1, 3} * U(:, [-3; -2] + 4 * leLIDs);
    vl = ul(:, 2 : 2 : end) ./ ul(:, 1 : 2 : end);

    Z12(:, leLIDs) = Z12(:, leLIDs) - J .* (bs.phitw_face{1, 3} * vl);

    % Take care of mass mtrix
    Z11 = (IME * Z11) ./ msh.elemJac(:, msh.LElems);
    Z12 = (IME * Z12) ./ msh.elemJac(:, msh.LElems);
    Z21 = (IME * Z21) ./ msh.elemJac(:, msh.LElems);
    Z22 = (IME * Z22) ./ msh.elemJac(:, msh.LElems);

    % viscous stress tensor
    TAU11 = tc.mu * (4 / 3 * Z11 - 2 / 3 * Z22);
    TAU12 = tc.mu * (Z12 +  Z21);
    TAU22 = tc.mu * (-2 / 3 * Z11 + 4 / 3 * Z22);
end

%**************************************************************************
%                       element contributions   
%**************************************************************************
if isViscous
    tau11 = bs.phi{1} * TAU11;
    tau12 = bs.phi{1} * TAU12;
    tau22 = bs.phi{1} * TAU22;
else
    tau11 = [];
    tau12 = [];
    tau22 = [];
end

% block-(1,1)
e11 = 0;
for n = 1 : 16
    if (n == 1 || n == 4)
        continue;
    end
    i = ceil(n / 4); j = n - (i - 1) * 4;
    ent1 = getJacobianEntry1(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), tau11, tau12, i, j, isViscous, tc);
    ent2 = getJacobianEntry2(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), tau12, tau22, i, j, isViscous, tc);
    
    s11 = e11 + 1; l = msh.nLElems * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s11 : e11) = -bs.phitw{2} * (repelem(msh.elemJxix( :, msh.LElems) .* ent1, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])) + ...
                     -bs.phitw{3} * (repelem(msh.elemJetay(:, msh.LElems) .* ent2, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
end

if isViscous
    % block-(1,2)
    i = 4;
    for j = 1 : 4
        ent1 = getJacobianEntry1(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), tau11, tau12, i, j + 4, isViscous, tc);
        ent2 = getJacobianEntry2(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), tau12, tau22, i, j + 4, isViscous, tc);

        s12 = e12 + 1; l = msh.nLElems * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        v12(s12 : e12) = -bs.phitw{2} * (repelem(msh.elemJxix( :, msh.LElems) .* ent1, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])) + ...
                         -bs.phitw{3} * (repelem(msh.elemJetay(:, msh.LElems) .* ent2, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    end

    % block-(2,1)
    e21 = 0;
    for n = 1 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i1 = 2 * i - 1; i2 = 2 * i; 
        if (i <= 2)
            [ent, flag] = getVelocityDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), i, j);
        elseif (i == 3)
            [ent, flag] = getScaledTemperatureDeriv(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), j, tc);
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
end

%**************************************************************************
%                        face contributions   
%**************************************************************************
% internal vertical face contributions
pl = getPressure(ul1(1 : bs.nfp, :), ul1(bs.nfp + 1 : 2 * bs.nfp, :), ul1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul1(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur1(1 : bs.nfp, :), ur1(bs.nfp + 1 : 2 * bs.nfp, :), ur1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur1(3 * bs.nfp + 1 : end, :), tc);
if isViscous
    tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs1);
    tau11r = bs.phi_face{1, 1} * TAU11(:, reLIDs1);
    tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs1);
    tau12r = bs.phi_face{1, 1} * TAU12(:, reLIDs1);    
else
    tau11l = [];
    tau11r = [];
    tau12l = [];
    tau12r = [];     
end

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getIntDissCoeAndDeriv(ul1(1 : bs.nfp, :), ur1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vl1(bs.nfp + 1 : 2 * bs.nfp, :), vr1(bs.nfp + 1 : 2 * bs.nfp, :), pl, pr, 1, tc);
    
% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    
    entl = getJacobianEntry1(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), tau11l, tau12l, i, j, isViscous, tc);
    entr = getJacobianEntry1(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), tau11r, tau12r, i, j, isViscous, tc);   
    dHl  = 0.5 * (entl - dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    dHr  = 0.5 * (entr - dalphar((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur1((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul1((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s11 = e11 + 1; l = 4 * nf1 * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;  
    v11(s11 : e11) = [bs.phitw_face{1, 2} * (repelem(J1 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                      bs.phitw_face{1, 2} * (repelem(J1 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(J1 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(J1 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];    
end

if isViscous
    % block-(1,2)
    i = 4; 
    for j = 1 : 4
        entl = getJacobianEntry1(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), tau11l, tau12l, i, j + 4, isViscous, tc);
        entr = getJacobianEntry1(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), tau11r, tau12r, i, j + 4, isViscous, tc);  
        dHl  = 0.5 * entl;
        dHr  = 0.5 * entr;

        s12 = e12 + 1; l = 4 * nf1 * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
        c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;  
        v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(J1 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                          bs.phitw_face{1, 2} * (repelem(J1 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                         -bs.phitw_face{1, 1} * (repelem(J1 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                         -bs.phitw_face{1, 1} * (repelem(J1 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];     
    end

    % block-(2,1)
    for n = 1 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i1 = 2 * i - 1; 
        if (i <= 2)
            entl         = getVelocityDeriv(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), i, j);
            [entr, flag] = getVelocityDeriv(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            entl         = getScaledTemperatureDeriv(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), j, tc);
            [entr, flag] = getScaledTemperatureDeriv(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), j, tc);
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
end

% internal horizontal face contributions
pl = getPressure(ul2(1 : bs.nfp, :), ul2(bs.nfp + 1 : 2 * bs.nfp, :), ul2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul2(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur2(1 : bs.nfp, :), ur2(bs.nfp + 1 : 2 * bs.nfp, :), ur2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur2(3 * bs.nfp + 1 : end, :), tc);
if isViscous
    tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs2);
    tau12r = bs.phi_face{1, 3} * TAU12(:, reLIDs2);
    tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs2);
    tau22r = bs.phi_face{1, 3} * TAU22(:, reLIDs2);   
else
    tau12l = [];
    tau12r = [];
    tau22l = [];
    tau22r = [];     
end

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getIntDissCoeAndDeriv(ul2(1 : bs.nfp, :), ur2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vl2(bs.nfp + 1 : 2 * bs.nfp, :), vr2(bs.nfp + 1 : 2 * bs.nfp, :), pl, pr, 2, tc);
              
% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    
    entl = getJacobianEntry2(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), tau12l, tau22l, i, j, isViscous, tc);
    entr = getJacobianEntry2(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), tau12r, tau22r, i, j, isViscous, tc);   
    dHl  = 0.5 * (entl - dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    dHr  = 0.5 * (entr - dalphar((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur2((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul2((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s11 = e11 + 1; l = 4 * nf2 * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;  
    v11(s11 : e11) = [bs.phitw_face{1, 4} * (repelem(J2 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                      bs.phitw_face{1, 4} * (repelem(J2 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), ...
                     -bs.phitw_face{1, 3} * (repelem(J2 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                     -bs.phitw_face{1, 3} * (repelem(J2 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2]))];    
end

if isViscous
    % block-(1,2)
    i = 4; 
    for j = 1 : 4
        entl = getJacobianEntry2(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), tau12l, tau22l, i, j + 4, isViscous, tc);
        entr = getJacobianEntry2(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), tau12r, tau22r, i, j + 4, isViscous, tc);  
        dHl  = 0.5 * entl;
        dHr  = 0.5 * entr;

        s12 = e12 + 1; l = 4 * nf2 * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + leLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1, (i - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;
        c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1, (j - 1) * msh.nLElems + leLIDs2 - 1, (j - 1) * msh.nLElems + reLIDs2 - 1] * bs.nb;  
        v12(s12 : e12) = [bs.phitw_face{1, 4} * (repelem(J2 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                          bs.phitw_face{1, 4} * (repelem(J2 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2])), ...
                         -bs.phitw_face{1, 3} * (repelem(J2 .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf2])), ...
                         -bs.phitw_face{1, 3} * (repelem(J2 .* dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf2]))];     
    end

    % block-(2,1)
    for n = 1 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i2 = 2 * i; 
        if (i <= 2)
            entl         = getVelocityDeriv(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), i, j);
            [entr, flag] = getVelocityDeriv(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            entl         = getScaledTemperatureDeriv(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), j, tc);
            [entr, flag] = getScaledTemperatureDeriv(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), j, tc);
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
end

%**************************************************************************
%                      boundary face contributions   
%**************************************************************************
% slip solid wall boundary face contributions
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);

% dissipation coefficient and its derivative 
[alpha, dalphal] = getBndDissCoeAndDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : 2 * bs.nfp, :), pl, 1, [-1; 0], tc);

% block-(1,1)
i = 2;
for j = 1 : 4
    entl = getJacobianEntry1(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), [], [], 2, j, false, tc);   
    dHl  = -entl + dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* ul(bs.nfp + 1 : 2 * bs.nfp, :) + alpha * (j == 2);
      
    s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 1} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]));    
end

if isViscous
    tau12l = bs.phi_face{1, 1} * TAU12(:, leLIDs);
    % block-(1,1)
    i = 4;
    for j = [1, 3]
        if (j == 1)
            dHl = tau12l .* vl(bs.nfp + 1 : 2 * bs.nfp, :) ./ ul(1 : bs.nfp, :);
        elseif (j == 3)
            dHl = -tau12l ./ ul(1 : bs.nfp, :);
        end
        
        s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
        r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
        v11(s11 : e11) = -bs.phitw_face{1, 1} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]));
    end
    
    % block-(1,2)
    i = 4;
    for j = 2 : 3
        dHl = -tc.mu * vl(bs.nfp + 1 : 2 * bs.nfp, :);
        
        s12 = e12 + 1; l = nf * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb;  
        v12(s12 : e12) = -bs.phitw_face{1, 1} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]));
    end  
    
    % block-(2,1)
    leIDs = leLIDs;
    for n = 5 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i1 = 2 * i - 1; 
        if (i == 2)
            [entl, flag] = getVelocityDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            [entl, flag] = getScaledTemperatureDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), j, tc);
        end
        if (flag == 0)
           continue; 
        end  

        s21 = e21 + 1; l = nf * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + ((i1 - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        v21(s21 : e21) = reshape(-2 * IME * bs.phitw_face{1, 1} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemLength(1, leIDs);
    end    
end

faceIDs = md.bndLFaces{3, 1};
leLIDs  = msh.faceElems(1, faceIDs); 
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);

% dissipation coefficient and its derivative 
[alpha, dalphal] = getBndDissCoeAndDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : 2 * bs.nfp, :), pl, 2, [0; -1], tc);

% block-(1,1)
i = 3;
for j = 1 : 4
    entl = getJacobianEntry2(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), [], [], 3, j, false, tc);   
    dHl  = -entl + dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* ul(2 * bs.nfp + 1 : 3 * bs.nfp, :) + alpha * (j == 3);
      
    s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 3} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf]));    
end

if isViscous
    % block-(1,1)
    tau12l = bs.phi_face{1, 3} * TAU12(:, leLIDs);
    i = 4;
    for j = 1 : 2
        if (j == 1)
            dHl = tau12l .* vl(1 : bs.nfp, :) ./ ul(1 : bs.nfp, :);
        elseif (j == 2)
            dHl = -tau12l ./ ul(1 : bs.nfp, :);
        end
        
        s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
        r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
        v11(s11 : e11) = -bs.phitw_face{1, 3} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf]));
    end    
    
    % block-(1,2)
    i = 4;
    for j = 2 : 3
        dHl = -tc.mu * vl(1 : bs.nfp, :);
        
        s12 = e12 + 1; l = nf * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb;  
        v12(s12 : e12) = -bs.phitw_face{1, 3} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf]));
    end
    
    % block-(2,1)
    leIDs = leLIDs;
    for n = [1 : 4, 9 : 12]
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i2 = 2 * i; 
        if (i == 1)
            [entl, flag] = getVelocityDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            [entl, flag] = getScaledTemperatureDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), j, tc);
        end
        if (flag == 0)
           continue; 
        end  

        s21 = e21 + 1; l = nf * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + ((i2 - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        v21(s21 : e21) = reshape(-2 * IME * bs.phitw_face{1, 3} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf])), [bs.nb2, nf]) ./ msh.elemLength(2, leIDs);
    end    
end

% supersonic outflow boundary face contributions
faceIDs = md.bndLFaces{2, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
if isViscous
    tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs);
    tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs);  
end

% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    entl = getJacobianEntry1(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tau11l, tau12l, i, j, isViscous, tc);   
    dHl  = entl;
      
    s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 2} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf]));    
end

if isViscous
    % block-(1,2)
    i = 4; 
    for j = 1 : 4
        entl = getJacobianEntry1(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tau11l, tau12l, i, j + 4, isViscous, tc);
        dHl  = entl;

        s12 = e12 + 1; l = nf * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb;  
        v12(s12 : e12) = bs.phitw_face{1, 2} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf]));     
    end  
    
    % block-(2,1)
    for n = 1 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i1 = 2 * i - 1; 
        if (i <= 2)
            [entl, flag] = getVelocityDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            [entl, flag] = getScaledTemperatureDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), j, tc);
        end
        if (flag == 0)
           continue; 
        end  

        s21 = e21 + 1; l = nf * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + ((i1 - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        v21(s21 : e21) = reshape(2 * IME * bs.phitw_face{1, 2} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemLength(1, leIDs);
    end    
end

faceIDs = md.bndLFaces{4, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
if isViscous
    tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs);  
    tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs);
end

% block-(1,1)
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    entl = getJacobianEntry2(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tau12l, tau22l, i, j, isViscous, tc);   
    dHl  = entl;
      
    s11 = e11 + 1; l = nf * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 4} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf]));    
end

if isViscous
    % block-(1,2)
    i = 4; 
    for j = 1 : 4
        entl = getJacobianEntry2(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tau12l, tau22l, i, j + 4, isViscous, tc); 
        dHl  = entl;

        s12 = e12 + 1; l = nf * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb;  
        v12(s12 : e12) = bs.phitw_face{1, 4} * (repelem(J .* dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf]));     
    end

    % block-(2,1)
    for n = 1 : 12
        i = ceil(n / 4); j = n - (i - 1) * 4;
        i2 = 2 * i; 
        if (i <= 2)
            [entl, flag] = getVelocityDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), i, j);
        elseif (i == 3)
            [entl, flag] = getScaledTemperatureDeriv(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), j, tc);
        end
        if (flag == 0)
           continue; 
        end  

        s21 = e21 + 1; l = nf * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + ((i2 - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        c21(s21 : e21) = ce1 + ((j  - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
        v21(s21 : e21) = reshape(2 * IME * bs.phitw_face{1, 4} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf])), [bs.nb2, nf]) ./ msh.elemLength(2, leIDs);
    end
end


v11(1 : e11) = dt * bt.diagCoe * v11(1 : e11); 
if isViscous
    v12(1 : e12) = dt * bt.diagCoe * v12(1 : e12); 
end
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
r11(e11 + 1 : end) = [];
c11(e11 + 1 : end) = [];
v11(e11 + 1 : end) = [];
v11 = mychop(v11, mtol);
M11 = sparse(r11, c11, v11, 4 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
if isViscous
    r12(e12 + 1 : end) = [];
    c12(e12 + 1 : end) = [];
    v12(e12 + 1 : end) = [];
    r21(e21 + 1 : end) = [];
    c21(e21 + 1 : end) = [];
    v21(e21 + 1 : end) = []; 
    v12 = mychop(v12, mtol);
    v21 = mychop(v21, mtol);

    M12 = sparse(r12, c12, v12, 4 * msh.nLElems * bs.nb, 6 * msh.nLElems * bs.nb);
    M21 = sparse(r21, c21, v21, 6 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
    Imp.Mat = M11 + M12 * M21;
else
    Imp.Mat = M11;
end

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
    p = getPressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);
        
    % Initialize row, column and value array for the KKT contributions
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);  
    rKKT = zeros(7 * msh.nLElems * bs.nb2 + 10 * msh.nLElems * bs.nb * bs.nep + 2 * KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;    
    
    % block (1,1)
    e = 0;
    for n = [1, 2, 3, 5, 6, 9, 11]
        i = ceil(n / 4); j = n - (i - 1) * 4;
        dp2 = getPressureHessian(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), i, j, tc);

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
        dp = getPressureDeriv(v(:, 1 : end / 2), v(:, end / 2 + 1 : end), j, tc);
        
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 4 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp, 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 5 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    [~, set2] = getActiveSets(u(:, 1 : 4 : end), p, Lamb, dU, dLamb, rKKT(stmp : etmp), cKKT(stmp : etmp), vKKT(stmp : etmp), KKT, newton);    
      
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
% Compute dissipation coefficient and its derivatives at internal faces
function [alpha, dalphal, dalphar] = getIntDissCoeAndDeriv(rhol, rhor, v1l, v1r, v2l, v2r, pl, pr, dir, tc)

% number of face points on each face and number of faces
[nfp, nf_] = size(rhol);

% dissipation coefficient
if (dir == 1)
    alphal = getEigenMax(rhol, v1l, pl, tc);
    alphar = getEigenMax(rhor, v1r, pr, tc);
elseif (dir == 2)
    alphal = getEigenMax(rhol, v2l, pl, tc);
    alphar = getEigenMax(rhor, v2r, pr, tc);    
end
indl = alphal >= alphar;
indr = ~indl;
alpha = alphar;
alpha(indl) = alphal(indl);

% derivative of dissipation coefficient with respect to left and right state
if (dir == 1)
    vnl = v1l;
    vnr = v1r;
    ind = nfp + 1 : 2 * nfp;
elseif (dir == 2)
    vnl = v2l;
    vnr = v2r;
    ind = 2 * nfp + 1 : 3 * nfp;
end
sigl = sign(vnl); 
sigr = sign(vnr);  
sigl(sigl == 0) = 1; 
sigr(sigr == 0) = 1;

dalphal = repmat(0.5 * sqrt(abs(tc.gamma ./ rhol ./ pl)), [4, 1]) .* [0.5 * (tc.gamma - 1) * (v1l.^2 + v2l.^2) - pl ./ rhol; (1 - tc.gamma) * v1l; (1 - tc.gamma) * v2l; repmat(tc.gamma - 1, [nfp, nf_])];
dalphar = repmat(0.5 * sqrt(abs(tc.gamma ./ rhor ./ pr)), [4, 1]) .* [0.5 * (tc.gamma - 1) * (v1r.^2 + v2r.^2) - pr ./ rhor; (1 - tc.gamma) * v1r; (1 - tc.gamma) * v2r; repmat(tc.gamma - 1, [nfp, nf_])];

tmp = sigl ./ rhol;
dalphal(1 : nfp, :) = dalphal(1 : nfp, :) - tmp .* vnl;
dalphal(ind, :)     = dalphal(ind, :)     + tmp;

tmp = sigr ./ rhor;
dalphar(1 : nfp, :) = dalphar(1 : nfp, :) - tmp .* vnr;
dalphar(ind, :)     = dalphar(ind, :)     + tmp;

dalphal(repmat(indr, [4, 1])) = 0;
dalphar(repmat(indl, [4, 1])) = 0;

end

% Compute dissipation coefficient and its derivatives at slip solid wall 
% boundary faces
function [alpha, dalphal] = getBndDissCoeAndDeriv(rhol, v1l, v2l, pl, dir, n, tc)

% number of face points on each face and number of faces
[nfp, nf_] = size(rhol);

% dissipation coefficient
if (dir == 1)
    alpha = getEigenMax(rhol, v1l, pl, tc);
elseif (dir == 2)
    alpha = getEigenMax(rhol, v2l, pl, tc);
end

% derivative of dissipation coefficient with respect to left state
if (dir == 1)
    vnl = v1l .* n(1, :);
    ind = nfp + 1 : 2 * nfp;
elseif (dir == 2)
    vnl = v2l .* n(2, :);
    ind = 2 * nfp + 1 : 3 * nfp;
end
sig = sign(vnl); 
sig(sig == 0) = 1; 

dalphal = repmat(0.5 * sqrt(abs(tc.gamma ./ rhol ./ pl)), [4, 1]) .* [0.5 * (tc.gamma - 1) * (v1l.^2 + v2l.^2) - pl ./ rhol; (1 - tc.gamma) * v1l; (1 - tc.gamma) * v2l; repmat(tc.gamma - 1, [nfp, nf_])];

tmp = sig ./ rhol;  
dalphal(1 : nfp, :) = dalphal(1 : nfp, :) - tmp .* vnl;
dalphal(ind, :)     = dalphal(ind, :)     + tmp .* n(dir, :);

end
