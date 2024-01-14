function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, isViscous, pr12, pc12, pv12, IME, tc, bs, bt, KKT, newton, mtol)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

% Evaluate physical quantities 
u = bs.phi{1} * U;
v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);

faceIDs1 = md.intLFaces{1, 3};
leIDs1   = msh.faceElems(1, faceIDs1);
reIDs1   = msh.faceElems(2, faceIDs1);
leLIDs1  = leIDs1;
reLIDs1  = reIDs1;
nf1      = length(faceIDs1);

ul1 = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs1), [3, nf1]);
ur1 = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs1), [3, nf1]);
vl1 = ul1(2, :) ./ ul1(1, :);
vr1 = ur1(2, :) ./ ur1(1, :);

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);

nBndLFaces = sum(md.nBndLFaces);
r11 = zeros((9 * (msh.nLElems + 4 * md.nIntLFaces + nBndLFaces) + 3 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

if isViscous
    r12 = zeros(length(pr12) + (msh.nLElems + 4 * md.nIntLFaces + nBndLFaces) * bs.nb2, 1);
    c12 = r12;
    v12 = r12;

    r21 = zeros(5 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
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
else
    tau  = [];
    taul = [];
    taur = [];
end
%**************************************************************************
%                         element contributions   
%**************************************************************************
% block-(1,1)
e11 = 0;
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;
    [ent, flag] = getJacobianEntry(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), tau, i, j, isViscous, tc);
    if (flag == 0)
       continue; 
    end
    
    s11 = e11 + 1; l = msh.nLElems * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s11 : e11) = -bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
end

if isViscous
    % block-(1,2)
    i = 3; j = 1;
    [ent, flag] = getJacobianEntry(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), tau, i, j + 3, isViscous, tc);

    s12 = e12 + 1; l = msh.nLElems * bs.nb2; e12 = e12 + l;
    r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v12(s12 : e12) = -bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));

    % block-(2,1)
    e21 = 0;
    for n = 1 : 6
        i = ceil(n / 3); j = n - (i - 1) * 3;
        if (i == 1)
            [ent, flag] = getVelocityDeriv(u(:, 1 : 3 : end), v, j);
        elseif (i == 2)
            [ent, flag] = getScaledTemperatureDeriv(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), j, tc);
        end
        if (flag == 0)
           continue; 
        end    

        s21 = e21 + 1; l = msh.nLElems * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        c21(s21 : e21) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        v21(s21 : e21) = reshape(-IME * bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems); 
    end
end

%**************************************************************************
%                       internal face contributions   
%**************************************************************************
pl = getPressure(ul1(1, :), ul1(2, :), ul1(3, :), tc);
pr = getPressure(ur1(1, :), ur1(2, :), ur1(3, :), tc);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getIntDissCoeAndDeriv(ul1(1, :), ur1(1, :), vl1, vr1, pl, pr, tc);
           
% block-(1,1)
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;
    entl = getJacobianEntry(ul1(1, :), vl1, ul1(3, :), taul, i, j, isViscous, tc);
    entr = getJacobianEntry(ur1(1, :), vr1, ur1(3, :), taur, i, j, isViscous, tc);
    dHl  = 0.5 * (entl - dalphal(j, :) .* (ur1(i, :) - ul1(i, :)) + alpha * (i == j));
    dHr  = 0.5 * (entr - dalphar(j, :) .* (ur1(i, :) - ul1(i, :)) - alpha * (i == j));

    s11 = e11 + 1; l = 4 * nf1 * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb; 
    v11(s11 : e11) = [bs.phitw_face{1, 2} * (repelem(dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                      bs.phitw_face{1, 2} * (repelem(dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];    
end

if isViscous
    % block-(1,2)
    i = 3; j = 1;
    entl = getJacobianEntry(ul1(1, :), vl1, ul1(3, :), taul, i, j + 3, isViscous, tc);
    entr = getJacobianEntry(ur1(1, :), vr1, ur1(3, :), taur, i, j + 3, isViscous, tc);
    dHl  = 0.5 * entl;
    dHr  = 0.5 * entr;

    s12 = e12 + 1; l = 4 * nf1 * bs.nb2; e12 = e12 + l;
    r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
    c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb; 
    v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                      bs.phitw_face{1, 2} * (repelem(dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(dHl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), ...
                     -bs.phitw_face{1, 1} * (repelem(dHr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1]))];

    % block-(2,1)
    for n = 1 : 6
        i = ceil(n / 3); j = n - (i - 1) * 3;
        if (i == 1)
            entl         = getVelocityDeriv(ul1(1, :), vl1, j);
            [entr, flag] = getVelocityDeriv(ur1(1, :), vr1, j);
        elseif (i == 2)
            entl         = getScaledTemperatureDeriv(ul1(1, :), vl1, ul1(3, :), j, tc);
            [entr, flag] = getScaledTemperatureDeriv(ur1(1, :), vr1, ur1(3, :), j, tc);
        end
        if (flag == 0)
           continue; 
        end  

        s21 = e21 + 1; l = 4 * nf1 * bs.nb2; e21 = e21 + l;    
        r21(s21 : e21) = re1 + [(i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + leLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1, (i - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
        c21(s21 : e21) = ce1 + [(j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1, (j - 1) * msh.nLElems + leLIDs1 - 1, (j - 1) * msh.nLElems + reLIDs1 - 1] * bs.nb;
        v21(s21 : e21) = [reshape( IME * bs.phitw_face{1, 2} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(:, leIDs1), ...
                          reshape( IME * bs.phitw_face{1, 2} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(:, leIDs1), ...
                          reshape(-IME * bs.phitw_face{1, 1} * (repelem(entl, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(:, reIDs1), ...
                          reshape(-IME * bs.phitw_face{1, 1} * (repelem(entr, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf1])), [bs.nb2, nf1]) ./ msh.elemLength(:, reIDs1)];
    end
end

%**************************************************************************
%                      boundary face contributions   
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

% dissipation coefficient and its derivative 
[alpha, dalphal] = getBndDissCoeAndDeriv(ul(1, :), ur(1, :), vl, vr, pl, pr, -1, 2, tc);

% block-(1,1)
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;  
    entl = getJacobianEntry(ul(1, :), vl, ul(3, :), [], i, j, false, tc);
    dHl  = 0.5 * (-entl - dalphal(j, :) .* (ur(i, :) - ul(i, :)) + alpha * (i == j));
    
    s11 = e11 + 1; l = bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 1} * (dHl .* bs.phi_face{1, 1});    
end

if isViscous
    % block-(1,2)
    i = 3; j = 1;
    dHl = -4 / 3 * tc.mu * tc.uL;

    s12 = e12 + 1; l = bs.nb2; e12 = e12 + l;
    r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v12(s12 : e12) = -bs.phitw_face{1, 1} * (dHl .* bs.phi_face{1, 1});    
end

% right boundary face contributions   
faceIDs = md.bndLFaces{2, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoR; tc.mR; tc.ER];
vl = ul(2, :) ./ ul(1, :);
vr = tc.uR;
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pR;

% dissipation coefficient and its derivative 
[alpha, dalphal] = getBndDissCoeAndDeriv(ul(1, :), ur(1, :), vl, vr, pl, pr, 1, 2, tc);

% block-(1,1)
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;
    
    entl = getJacobianEntry(ul(1, :), vl, ul(3, :), [], i, j, false, tc);
    dHl  = 0.5 * (entl - dalphal(j, :) .* (ur(i, :) - ul(i, :)) + alpha * (i == j));
    
    s11 = e11 + 1; l = bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v11(s11 : e11) = bs.phitw_face{1, 2} * (dHl .* bs.phi_face{1, 2});    
end

if isViscous
    % block-(1,2)
    i = 3; j = 1;
    dHl = -4 / 3 * tc.mu * tc.uR;

    s12 = e12 + 1; l = bs.nb2; e12 = e12 + l;
    r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems + leLIDs - 1) * bs.nb;
    c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems + leLIDs - 1) * bs.nb; 
    v12(s12 : e12) = bs.phitw_face{1, 2} * (dHl .* bs.phi_face{1, 2});  
end
                    
v11(1 : e11) = dt * bt.diagCoe * v11(1 : e11); 
if isViscous
    v12(1 : e12) = dt * bt.diagCoe * v12(1 : e12); 
end
%**************************************************************************
%                        mass matrix contributions   
%**************************************************************************
s11 = e11 + 1; l = 3 * msh.nLElems * bs.nb2; e11 = e11 + l;
r11(s11 : e11) = re1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
c11(s11 : e11) = ce1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
v11(s11 : e11) = repmat(msh.elemJac(:, msh.LElems) .* reshape(bs.phitw{1} * bs.phi{1}, [bs.nb2, 1]), [1, 3]);

%**************************************************************************
%                   global matrix for DG discretization
%**************************************************************************
% Assemble
r11(e11 + 1 : end) = [];
c11(e11 + 1 : end) = [];
v11(e11 + 1 : end) = [];
v11 = mychop(v11, mtol);
M11 = sparse(r11, c11, v11, 3 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb);
if isViscous
    r12(e12 + 1 : end) = [];
    c12(e12 + 1 : end) = [];
    v12(e12 + 1 : end) = [];
    r21(e21 + 1 : end) = [];
    c21(e21 + 1 : end) = [];
    v21(e21 + 1 : end) = [];    
    v12 = mychop(v12, mtol);
    v21 = mychop(v21, mtol);

    M12 = sparse(r12, c12, v12, 3 * msh.nLElems * bs.nb, 2 * msh.nLElems * bs.nb);
    M21 = sparse(r21, c21, v21, 2 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb);
    Imp.Mat = M11 + M12 * M21;
else
    Imp.Mat = M11;
end

%**************************************************************************
%                       KKT constraints contribution 
%**************************************************************************
% Add positivity constraints contribution
if (KKT.positivity)      
    if (bt.form == 1)
        cdt  = 1;
        cdt2 = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;
        cdt2 = cdt * cdt;
    end  
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
        
    % Initialize row, column and value array for the KKT contributions
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);  
    rKKT = zeros(4 * msh.nLElems * bs.nb2 + 8 * msh.nLElems * bs.nb * bs.nep + 2 * KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;    
    
    % block (1,1)
    e = 0;
    for n = [1, 2, 4, 5]
        i = ceil(n / 3); j = n - (i - 1) * 3;
        dp2 = getPressureHessian(u(:, 1 : 3 : end), v, i, j, tc);

        s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
        rKKT(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        cKKT(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
    end    
    
    % KKT-block (3,1)    
    s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
    stmp = s;
    rKKT(s : e) = 3 * KKT.NDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
    cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
    vKKT(s : e) = repmat(cdt * bs.phi{1}(:), [msh.nLElems, 1]);
    for j = 1 : 3
        dp = getPressureDeriv(v, j, tc);
        
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 3 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp, 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 4 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
    [~, set2] = getActiveSets(u(:, 1 : 3 : end), p, Lamb, dU, dLamb, rKKT(stmp : etmp), cKKT(stmp : etmp), vKKT(stmp : etmp), KKT, newton);    
      
    % KKT-block (3,3)
    s = e + 1; l = 2 * KKT.LDOF; e = e + l;
    rKKT(s : e) = 3 * KKT.NDOF + (1 : 2 * KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    [rDG, cDG, vDG] = find(Imp.Mat);
    vKKT = mychop(vKKT, mtol);
    Imp.Mat = sparse([rDG; rKKT], [cDG; cKKT], [vDG; vKKT], 3 * KKT.NDOF + 2 * KKT.LDOF, 3 * KKT.NDOF + 2 * KKT.LDOF);
    Imp.Mat([false(3 * KKT.NDOF, 1); set2], 1 : 3 * KKT.NDOF) = 0;
     
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix
        Imp.GT  = Imp.Mat';
        Imp.Mat = Imp.GT * Imp.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(Imp.Mat, 1));
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        Imp.GT = speye(size(Imp.Mat, 1));
    end
else
    Imp.GT  = speye(3 * KKT.NDOF);
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
function [alpha, dalphal, dalphar] = getIntDissCoeAndDeriv(rhol, rhor, vl, vr, pl, pr, tc)
    
% number of faces 
nf = size(rhol, 2);

% dissipation coefficient
alphal = getEigenMax(rhol, vl, pl, tc);
alphar = getEigenMax(rhor, vr, pr, tc);
alpha  = max(alphal, alphar);

ind1 = alphal >= alphar;
ind2 = ~ind1;
nf1  = sum(ind1);
nf2  = nf - nf1;

% derivative of dissipation coefficient with respect to left and right state
dalphal = zeros(3, nf);
dalphar = zeros(3, nf);

dalphal(:, ind1) = 0.5 * sqrt(abs(tc.gamma ./ rhol(ind1) ./ pl(ind1))) .* [0.5 * (tc.gamma - 1) * vl(ind1).^2 - pl(ind1) ./ rhol(ind1); (1 - tc.gamma) * vl(ind1); repmat(tc.gamma - 1, [1, nf1])];
dalphar(:, ind2) = 0.5 * sqrt(abs(tc.gamma ./ rhor(ind2) ./ pr(ind2))) .* [0.5 * (tc.gamma - 1) * vr(ind2).^2 - pr(ind2) ./ rhor(ind2); (1 - tc.gamma) * vr(ind2); repmat(tc.gamma - 1, [1, nf2])];

sig = sign(vl(ind1)); sig(sig == 0) = 1;
dalphal(1 : 2, ind1) = dalphal(1 : 2, ind1) + sig ./ rhol(ind1) .* [-vl(ind1); ones(1, nf1)];

sig = sign(vr(ind2)); sig(sig == 0) = 1;
dalphar(1 : 2, ind2) = dalphar(1 : 2, ind2) + sig ./ rhor(ind2) .* [-vr(ind2); ones(1, nf2)];

end

% Compute dissipation coefficient and its derivatives at boundary faces
function [alpha, dalphal] = getBndDissCoeAndDeriv(rhol, rhor, vl, vr, pl, pr, n, bndType, tc)

% number of faces 
nf = size(rhol, 2);

% dissipation coefficient
alphal = getEigenMax(rhol, vl, pl, tc);
alphar = getEigenMax(rhor, vr, pr, tc);
alpha  = max(alphal, alphar);

ind1 = alphal >= alphar;
ind2 = ~ind1;
nf1  = sum(ind1);
nf2  = nf - nf1;

% derivative of dissipation coefficient with respect to left state
dalphal = zeros(3, nf);

if (nf1 > 0)
    dalphal(:, ind1) = 0.5 * sqrt(abs(tc.gamma ./ rhol(ind1) ./ pl(ind1))) .* [0.5 * (tc.gamma - 1) * vl(ind1).^2 - pl(ind1) ./ rhol(ind1); (1 - tc.gamma) * vl(ind1); repmat(tc.gamma - 1, [1, nf1])];

    vnl = vl(ind1) * n(ind1);
    sig = sign(vnl); sig(sig == 0) = 1;
    dalphal(1 : 2, ind1) = dalphal(1 : 2, ind1) + sig ./ rhol(ind1) .* [-vnl; n(ind1)];
end

if (nf2 > 0)
    switch bndType
        % Dirichlet boundary
        case 2
            dalphal(:, ind2) = 0;
    end
end
end
