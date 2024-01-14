% Compute the DG residual
function res = computeDGResidual(msh, md, U, isViscous, IME, tc, bs)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(U, [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

% Evaluate physical quantities 
u = bs.phi{1} * U; 
v = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];

% vertical faces
faceIDs = md.intLFaces{1, 3};
leLIDs1 = msh.faceElems(1, faceIDs);
reLIDs1 = msh.faceElems(2, faceIDs);
J1      = msh.faceJac(:, faceIDs);
nf1     = length(faceIDs);

ul1 = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs1), [4 * bs.nfp, nf1]);
ur1 = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs1), [4 * bs.nfp, nf1]);
vl1 = ul1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul1(1 : bs.nfp, :), [2, 1]);
vr1 = ur1(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur1(1 : bs.nfp, :), [2, 1]);

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

%**************************************************************************
%                     Compute viscous stress tensor
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
%                           Compute heat flux
%************************************************************************** 
if isViscous
    % element contributions
    sT = getScaledTemperature(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), tc);
    Q1 = -msh.elemJxix( :, msh.LElems) .* (bs.phitw{2} * sT);
    Q2 = -msh.elemJetay(:, msh.LElems) .* (bs.phitw{3} * sT);

    % internal vertical face contributions
    sTl = getScaledTemperature(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), vl1(bs.nfp + 1 : end, :), ul1(3 * bs.nfp + 1 : end, :), tc);
    sTr = getScaledTemperature(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), vr1(bs.nfp + 1 : end, :), ur1(3 * bs.nfp + 1 : end, :), tc);

    sT_hat = 0.5 * (sTl + sTr);
    Q1(:, leLIDs1) = Q1(:, leLIDs1) + J1 .* (bs.phitw_face{1, 2} * sT_hat);
    Q1(:, reLIDs1) = Q1(:, reLIDs1) - J1 .* (bs.phitw_face{1, 1} * sT_hat);

    % internal horizontal face contributions
    sTl = getScaledTemperature(ul2(1 : bs.nfp, :), vl2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), ul2(3 * bs.nfp + 1 : end, :), tc);
    sTr = getScaledTemperature(ur2(1 : bs.nfp, :), vr2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), ur2(3 * bs.nfp + 1 : end, :), tc);

    sT_hat = 0.5 * (sTl + sTr);
    Q2(:, leLIDs2) = Q2(:, leLIDs2) + J2 .* (bs.phitw_face{1, 4} * sT_hat);
    Q2(:, reLIDs2) = Q2(:, reLIDs2) - J2 .* (bs.phitw_face{1, 3} * sT_hat);

    % slip solid wall boundary face contributions
    faceIDs = md.bndLFaces{1, 1};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul  = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
    vl  = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
    sTl = getScaledTemperature(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tc);

    Q1(:, leLIDs) = Q1(:, leLIDs) - J .* (bs.phitw_face{1, 1} * sTl);

    faceIDs = md.bndLFaces{3, 1};
    leLIDs  = msh.faceElems(1, faceIDs); 
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul  = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
    vl  = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
    sTl = getScaledTemperature(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tc);

    Q2(:, leLIDs) = Q2(:, leLIDs) - J .* (bs.phitw_face{1, 3} * sTl);

    % supersonic outflow boundary face contributions
    faceIDs = md.bndLFaces{2, 2};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul  = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
    vl  = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
    sTl = getScaledTemperature(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tc);

    Q1(:, leLIDs) = Q1(:, leLIDs) + J .* (bs.phitw_face{1, 2} * sTl);

    faceIDs = md.bndLFaces{4, 2};
    leLIDs  = msh.faceElems(1, faceIDs);
    J       = msh.faceJac(:, faceIDs);
    nf      = length(faceIDs);

    ul  = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
    vl  = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
    sTl = getScaledTemperature(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), tc);

    Q2(:, leLIDs) = Q2(:, leLIDs) + J .* (bs.phitw_face{1, 4} * sTl);

    % Take care of mass mtrix
    Q1 = (IME * Q1) ./ msh.elemJac(:, msh.LElems);
    Q2 = (IME * Q2) ./ msh.elemJac(:, msh.LElems);
end

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
if (bs.type == 100)
    res = zeros(1, 4 * msh.nLElems);
else
    p = getPressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);
    F = getF(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);
    G = getG(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);
    if isViscous
        tau11 = bs.phi{1} * TAU11;
        tau12 = bs.phi{1} * TAU12;
        tau22 = bs.phi{1} * TAU22;
        q1    = bs.phi{1} * Q1;
        q2    = bs.phi{1} * Q2;
    
        F(bs.nep + 1 : 2 * bs.nep, :)     = F(bs.nep + 1 : 2 * bs.nep, :)     - tau11;
        F(2 * bs.nep + 1 : 3 * bs.nep, :) = F(2 * bs.nep + 1 : 3 * bs.nep, :) - tau12;
        F(3 * bs.nep + 1 : end, :)        = F(3 * bs.nep + 1 : end, :)        - tau11 .* v(:, 1 : end / 2) - tau12 .* v(:, end / 2 + 1 : end) + q1;
        G(bs.nep + 1 : 2 * bs.nep, :)     = G(bs.nep + 1 : 2 * bs.nep, :)     - tau12;
        G(2 * bs.nep + 1 : 3 * bs.nep, :) = G(2 * bs.nep + 1 : 3 * bs.nep, :) - tau22;
        G(3 * bs.nep + 1 : end, :)        = G(3 * bs.nep + 1 : end, :)        - tau12 .* v(:, 1 : end / 2) - tau22 .* v(:, end / 2 + 1 : end) + q2;
    end
    
    res = -bs.phitw{2} * reshape(msh.elemJxix( :, msh.LElems) .* F, [bs.nep, 4 * msh.nLElems]) ...                  
         -bs.phitw{3} * reshape(msh.elemJetay(:, msh.LElems) .* G, [bs.nep, 4 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
% internal vertical face contributions
pl = getPressure(ul1(1 : bs.nfp, :), ul1(bs.nfp + 1 : 2 * bs.nfp, :), ul1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul1(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur1(1 : bs.nfp, :), ur1(bs.nfp + 1 : 2 * bs.nfp, :), ur1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur1(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul1(1 : bs.nfp, :), ul1(bs.nfp + 1 : 2 * bs.nfp, :), ul1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul1(3 * bs.nfp + 1 : end, :), pl); 
Fr = getF(ur1(1 : bs.nfp, :), ur1(bs.nfp + 1 : 2 * bs.nfp, :), ur1(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur1(3 * bs.nfp + 1 : end, :), pr);
if isViscous
    tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs1);
    tau11r = bs.phi_face{1, 1} * TAU11(:, reLIDs1);
    tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs1);
    tau12r = bs.phi_face{1, 1} * TAU12(:, reLIDs1);
    ql     = bs.phi_face{1, 2} * Q1(:, leLIDs1);
    qr     = bs.phi_face{1, 1} * Q1(:, reLIDs1);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11l;
    Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11r;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12l;
    Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12r;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau11l .* vl1(1 : bs.nfp, :) - tau12l .* vl1(bs.nfp + 1 : end, :) + ql;
    Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - tau11r .* vr1(1 : bs.nfp, :) - tau12r .* vr1(bs.nfp + 1 : end, :) + qr;
end

% dissipation coefficient 
alpha = max(getEigenMax(ul1(1 : bs.nfp, :), vl1(1 : bs.nfp, :), pl, tc), getEigenMax(ur1(1 : bs.nfp, :), vr1(1 : bs.nfp, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur1 - ul1));
res(:, (-3 : 0)' + 4 * leLIDs1) = res(:, (-3 : 0)' + 4 * leLIDs1) + bs.phitw_face{1, 2} * reshape(J1 .* F_hat, [bs.nfp, 4 * nf1]);
res(:, (-3 : 0)' + 4 * reLIDs1) = res(:, (-3 : 0)' + 4 * reLIDs1) - bs.phitw_face{1, 1} * reshape(J1 .* F_hat, [bs.nfp, 4 * nf1]);

% internal horizontal face contributions
pl = getPressure(ul2(1 : bs.nfp, :), ul2(bs.nfp + 1 : 2 * bs.nfp, :), ul2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul2(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur2(1 : bs.nfp, :), ur2(bs.nfp + 1 : 2 * bs.nfp, :), ur2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur2(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul2(1 : bs.nfp, :), ul2(bs.nfp + 1 : 2 * bs.nfp, :), ul2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul2(3 * bs.nfp + 1 : end, :), pl); 
Fr = getG(ur2(1 : bs.nfp, :), ur2(bs.nfp + 1 : 2 * bs.nfp, :), ur2(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur2(3 * bs.nfp + 1 : end, :), pr);
if isViscous
    tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs2);
    tau12r = bs.phi_face{1, 3} * TAU12(:, reLIDs2);
    tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs2);
    tau22r = bs.phi_face{1, 3} * TAU22(:, reLIDs2);
    ql     = bs.phi_face{1, 4} * Q2(:, leLIDs2);
    qr     = bs.phi_face{1, 3} * Q2(:, reLIDs2);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12l;
    Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12r;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22l;
    Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22r;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau12l .* vl2(1 : bs.nfp, :) - tau22l .* vl2(bs.nfp + 1 : end, :) + ql;
    Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - tau12r .* vr2(1 : bs.nfp, :) - tau22r .* vr2(bs.nfp + 1 : end, :) + qr;
end

% dissipation coefficient 
alpha = max(getEigenMax(ul2(1 : bs.nfp, :), vl2(bs.nfp + 1 : end, :), pl, tc), getEigenMax(ur2(1 : bs.nfp, :), vr2(bs.nfp + 1 : end, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur2 - ul2));
res(:, (-3 : 0)' + 4 * leLIDs2) = res(:, (-3 : 0)' + 4 * leLIDs2) + bs.phitw_face{1, 4} * reshape(J2 .* F_hat, [bs.nfp, 4 * nf2]);
res(:, (-3 : 0)' + 4 * reLIDs2) = res(:, (-3 : 0)' + 4 * reLIDs2) - bs.phitw_face{1, 3} * reshape(J2 .* F_hat, [bs.nfp, 4 * nf2]);

%**************************************************************************   
%               boundary face contributions for DG residual
%************************************************************************** 
% slip solid wall boundary face contributions
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
ur = ul;
ur(bs.nfp + 1 : 2 * bs.nfp, :) = -ul(bs.nfp + 1 : 2 * bs.nfp, :);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getF(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pl);
if isViscous
    tau11l = bs.phi_face{1, 1} * TAU11(:, leLIDs);
    tau12l = bs.phi_face{1, 1} * TAU12(:, leLIDs);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11l;
    Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11l;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12l;
    Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12l;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau12l .* vl(bs.nfp + 1 : end, :);
    Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - tau12l .* vl(bs.nfp + 1 : end, :);
end

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), pl, tc);

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{3, 1};
leLIDs  = msh.faceElems(1, faceIDs); 
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
ur = ul;
ur(2 * bs.nfp + 1 : 3 * bs.nfp, :) = -ul(2 * bs.nfp + 1 : 3 * bs.nfp, :);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getG(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pl);
if isViscous
    tau12l = bs.phi_face{1, 3} * TAU12(:, leLIDs);
    tau22l = bs.phi_face{1, 3} * TAU22(:, leLIDs);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12l;
    Fr(bs.nfp + 1 : 2 * bs.nfp, :)     = Fr(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12l;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22l;
    Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fr(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22l;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau12l .* vl(1 : bs.nfp, :);
    Fr(3 * bs.nfp + 1 : end, :)        = Fr(3 * bs.nfp + 1 : end, :)        - tau12l .* vl(1 : bs.nfp, :);
end

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), pl, tc);

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% supersonic outflow boundary face contributions
faceIDs = md.bndLFaces{2, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
if isViscous
    tau11l = bs.phi_face{1, 2} * TAU11(:, leLIDs);
    tau12l = bs.phi_face{1, 2} * TAU12(:, leLIDs);
    ql     = bs.phi_face{1, 2} * Q1(:, leLIDs);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau11l;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau12l;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau11l .* vl(1 : bs.nfp, :) - tau12l .* vl(bs.nfp + 1 : end, :) + ql;
end

F_hat = Fl;
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{4, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
if isViscous
    tau12l = bs.phi_face{1, 4} * TAU12(:, leLIDs);
    tau22l = bs.phi_face{1, 4} * TAU22(:, leLIDs);
    ql     = bs.phi_face{1, 4} * Q2(:, leLIDs);

    Fl(bs.nfp + 1 : 2 * bs.nfp, :)     = Fl(bs.nfp + 1 : 2 * bs.nfp, :)     - tau12l;
    Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) = Fl(2 * bs.nfp + 1 : 3 * bs.nfp, :) - tau22l;
    Fl(3 * bs.nfp + 1 : end, :)        = Fl(3 * bs.nfp + 1 : end, :)        - tau12l .* vl(1 : bs.nfp, :) - tau22l .* vl(bs.nfp + 1 : end, :) + ql;
end

F_hat = Fl;
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

%**************************************************************************   
%                            DG residual 
%************************************************************************** 
res = res(:, (1 : 4 : end)' + (0 : 3));
res = res(:);

end
