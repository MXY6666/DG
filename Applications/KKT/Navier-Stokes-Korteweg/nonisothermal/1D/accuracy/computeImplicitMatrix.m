function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, A, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, bs, bt, inviscidFlux, KKT, newton, mtol)

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

faceIDs = md.intLFaces{1, 3};
leIDs   = msh.faceElems(1, faceIDs);
reIDs   = msh.faceElems(2, faceIDs);
leLIDs  = leIDs;
reLIDs  = reIDs;
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
vl = ul(2, :) ./ ul(1, :);
vr = ur(2, :) ./ ur(1, :);

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);

r11 = zeros(((msh.nLElems + 4 * md.nIntLFaces) * 9 + 3 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

r12 = zeros(length(pr12) + 6 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c12 = r12;
v12 = r12;

r21 = zeros(length(pr21) + 5 * (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c21 = r21;
v21 = r21;

r22i = zeros(length(pr22i) + (msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
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
%                       element contributions   
%**************************************************************************
w  = bs.phi{1} * W;
r  = bs.phi{1} * R;
T  = getTemperature(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, tc);
p  = getPressure(u(:, 1 : 3 : end), T);
xi = getKortewegStressTensor(u(:, 1 : 3 : end), w, r, tc);
dT = getTemperatureDeriv(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, tc);
dp = getPressureDeriv(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, tc);

% block-(1,1)
e11 = 0;
for n = 1 : 9
    if (n == 1 || n == 3)
        continue;
    end    
    i = ceil(n / 3); j = n - (i - 1) * 3;
    ent = getJacobianEntry(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), p, w, r, z, tau, xi, dp, i, j, tc);
    
    s11 = e11 + 1; l = msh.nLElems * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s11 : e11) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s11 : e11) = -bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));     
end

% block-(1,2)
for n = 1 : 12   
    i = ceil(n / 4); j = n - (i - 1) * 4;
    [ent, flag] = getJacobianEntry(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), p, w, r, z, tau, xi, dp, i, j + 3, tc);
    
    if (flag == 2)
        s12 = e12 + 1; l = msh.nLElems * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        c12(s12 : e12) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        v12(s12 : e12) = -bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
    end    
end

% block-(2,1)
for n = 1 : 6
    i = ceil(n / 3); j = n - (i - 1) * 3; i = i + 2;
    if (i == 3)
        [ent, flag] = getVelocityDeriv(u(:, 1 : 3 : end), v, j);
    elseif (i == 4)
        ent = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dT(:, :, j);
        flag = 2;
    end
    if (flag == 0)
       continue; 
    end    
    
    s21 = e21 + 1; l = msh.nLElems * bs.nb2; e21 = e21 + l;    
    r21(s21 : e21) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c21(s21 : e21) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v21(s21 : e21) = reshape(-IME * bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems); 
end

% inverse of block-(2,2)
ent = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dT(:, :, 4);

s22 = e22 + 1; l = msh.nLElems * bs.nb2; e22 = e22 + l;    
r22i(s22 : e22) = re1 + (3 * msh.nLElems : 4 * msh.nLElems - 1) * bs.nb;
c22i(s22 : e22) = ce1 + (0 : msh.nLElems - 1) * bs.nb;
v22i(s22 : e22) = reshape(IME * bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems);

%**************************************************************************
%                        face contributions   
%**************************************************************************
wl  = bs.phi_face{1, 2} * W(:, leLIDs);
wr  = bs.phi_face{1, 1} * W(:, reLIDs);
rl  = bs.phi_face{1, 2} * R(:, leLIDs); 
rr  = bs.phi_face{1, 1} * R(:, reLIDs);
Tl  = getTemperature(ul(1, :), vl, ul(3, :), wl, tc);
Tr  = getTemperature(ur(1, :), vr, ur(3, :), wr, tc);
pl  = getPressure(ul(1, :), Tl);
pr  = getPressure(ur(1, :), Tr);
xil = getKortewegStressTensor(ul(1, :), wl, rl, tc);
xir = getKortewegStressTensor(ur(1, :), wr, rr, tc);
dTl = getTemperatureDeriv(ul(1, :), vl, ul(3, :), wl, tc);
dTr = getTemperatureDeriv(ur(1, :), vr, ur(3, :), wr, tc);
dpl = getPressureDeriv(ul(1, :), vl, ul(3, :), wl, tc);
dpr = getPressureDeriv(ur(1, :), vr, ur(3, :), wr, tc);

% dissipation coefficient and its derivatives 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul(1, :), ur(1, :), vl, vr, Tl, Tr, pl, pr, dTl, dTr, dpl, dpr, tc, inviscidFlux);
        
% block-(1,1)
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;
    entl = getJacobianEntry(ul(1, :), vl, ul(3, :), pl, wl, rl, zl, taul, xil, dpl, i, j, tc);
    entr = getJacobianEntry(ur(1, :), vr, ur(3, :), pr, wr, rr, zr, taur, xir, dpr, i, j, tc);
    H_ul = 0.5 * (entl - dalphal(j, :) .* (ur(i, :) - ul(i, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar(j, :) .* (ur(i, :) - ul(i, :)) - alpha * (i == j));
    
    s11 = e11 + 1; l = 4 * nf * bs.nb2; e11 = e11 + l;
    r11(s11 : e11) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    c11(s11 : e11) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb; 
    v11(s11 : e11) = [bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                      bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
                     -bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                     -bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];    
end

% block-(1,2)
for n = 1 : 12
    i = ceil(n / 4); j = n - (i - 1) * 4;
    entl         = getJacobianEntry(ul(1, :), vl, ul(3, :), pl, wl, rl, zl, taul, xil, dpl, i, j + 3, tc);
    [entr, flag] = getJacobianEntry(ur(1, :), vr, ur(3, :), pr, wr, rr, zr, taur, xir, dpr, i, j + 3, tc);
    
    if (j == 1)       
        H_ul = 0.5 * (entl - dalphal(4, :) .* (ur(i, :) - ul(i, :)));
        H_ur = 0.5 * (entr - dalphar(4, :) .* (ur(i, :) - ul(i, :)));
        
        s12 = e12 + 1; l = 4 * nf * bs.nb2; e12 = e12 + l;
        r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
        c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
        v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                          bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
                         -bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                         -bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];        
    else
        if (flag == 2)
            H_ul = 0.5 * entl;
            H_ur = 0.5 * entr;
            
            s12 = e12 + 1; l = 4 * nf * bs.nb2; e12 = e12 + l;
            r12(s12 : e12) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
            c12(s12 : e12) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb; 
            v12(s12 : e12) = [bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                              bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
                             -bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                             -bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];            
        end
    end
end

% block-(2,1)
for n = 1 : 6
    i = ceil(n / 3); j = n - (i - 1) * 3; i = i + 2;
    if (i == 3)
        entl         = getVelocityDeriv(ul(1, :), vl, j);
        [entr, flag] = getVelocityDeriv(ur(1, :), vr, j);
    elseif (i == 4)
        entl = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, j);
        entr = -8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, j);
        flag = 2;
    end
    if (flag == 0)
       continue; 
    end  
    H_ul = 0.5 * entl;
    H_ur = 0.5 * entr;
    
    s21 = e21 + 1; l = 4 * nf * bs.nb2; e21 = e21 + l;    
    r21(s21 : e21) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    c21(s21 : e21) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    v21(s21 : e21) = [reshape( IME * bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      reshape( IME * bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      reshape(-IME * bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                      reshape(-IME * bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)];
end

% inverse of block-(2,2)
H_ul = -0.5 * 8 * tc.Cv / (3 * tc.We * tc.Pr) * dTl(:, :, 4);
H_ur = -0.5 * 8 * tc.Cv / (3 * tc.We * tc.Pr) * dTr(:, :, 4);

s22 = e22 + 1; l = 4 * nf * bs.nb2; e22 = e22 + l;    
r22i(s22 : e22) = re1 + [3 * msh.nLElems + leLIDs - 1, 3 * msh.nLElems + leLIDs - 1, 3 * msh.nLElems + reLIDs - 1, 3 * msh.nLElems + reLIDs - 1] * bs.nb;
c22i(s22 : e22) = ce1 + [leLIDs - 1, reLIDs - 1, leLIDs - 1, reLIDs - 1] * bs.nb;
v22i(s22 : e22) = [reshape(-IME * bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                   reshape(-IME * bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                   reshape( IME * bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                   reshape( IME * bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)];
                    
v11(1 : e11) = dt * bt.diagCoe * v11(1 : e11); 
v12(1 : e12) = dt * bt.diagCoe * v12(1 : e12); 
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
M11  = sparse(r11,  c11,  v11,  3 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb);
M12  = sparse(r12,  c12,  v12,  3 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
M21  = sparse(r21,  c21,  v21,  4 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb);
M22i = sparse(r22i, c22i, v22i, 4 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);
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
    
    dp2 = getPressureHessian(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, 1, 4, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp11 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, 4, 4, tc);
    vtmp = -cdt2 * bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp12 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    Mtmp1 = Mtmp11 * A + A' * (Mtmp11 + Mtmp12 * A);
    [rKKT1, cKKT1, vKKT1] = find(Mtmp1);
    
    rtmp = re2 + (0 : msh.nLElems - 1) * bs.nep;
    ctmp = ce2 + (0 : msh.nLElems - 1) * bs.nb;  
    vtmp = repelem(dp(:, :, 4), 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    
    Mtmp2 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.LDOF, KKT.NDOF) * A;
    [rKKT2, cKKT2, vKKT2] = find(Mtmp2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    % Initialize row, column and value array for the KKT contributions
    rKKT = zeros(length(vKKT1) + 2 * length(vKKT2) + 6 * msh.nLElems * bs.nb2 + 8 * msh.nLElems * bs.nb * bs.nep + 2 * KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;    
    
    e = length(vKKT1) + 2 * length(vKKT2);
    rKKT(1 : e) = [rKKT1; 3 * KKT.NDOF + KKT.LDOF + rKKT2; cKKT2];
    cKKT(1 : e) = [cKKT1; cKKT2; 3 * KKT.NDOF + KKT.LDOF + rKKT2];
    vKKT(1 : e) = [vKKT1; vKKT2; -vKKT2];
    
    % block (1,1)
    for n = [1, 2, 3, 5, 6, 9]
        i = ceil(n / 4); j = n - (i - 1) * 4;
        dp2 = getPressureHessian(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), w, i, j, tc);

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
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 3 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp(:, :, j), 1, bs.nb) .* repmat(cdt * bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 4 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    [~, set2] = getActiveSets(u(:, 1 : 3 : end), p, Lamb, dU, dLamb, [rKKT(stmp : etmp); 3 * KKT.NDOF + KKT.LDOF + rKKT2], [cKKT(stmp : etmp); cKKT2], [vKKT(stmp : etmp); vKKT2], KKT, newton);    
      
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
% Compute dissipation coefficient and its derivatives
function [alpha, dalphal, dalphar] = getDissCoeAndDeriv(rhol, rhor, vl, vr, Tl, Tr, pl, pr, dTl, dTr, dpl, dpr, tc, inviscidFlux)
    
% number of faces
nf = size(rhol, 2);

% dissipation coefficient at each face
betal = 2 * (pl ./ rhol + 4 * Tl .* (3 + tc.Cv * (2 * rhol - 3)) ./ (tc.Cv * (rhol - 3).^2));
betar = 2 * (pr ./ rhor + 4 * Tr .* (3 + tc.Cv * (2 * rhor - 3)) ./ (tc.Cv * (rhor - 3).^2));

alphal = max(abs(vl - sqrt(betal)), abs(vl + sqrt(betal)));
alphar = max(abs(vr - sqrt(betar)), abs(vr + sqrt(betar)));
alpha  = max(alphal, alphar);

if (inviscidFlux == 1)
    % global dissipation coefficient
    [alpha, ind] = max(alpha);       
    
    % physical quantities evaluated at the point where global dissipation 
    % coefficient is get
    if (alphal(ind) >= alphar(ind))
        rho  = rhol(ind);
        v    = vl(ind); 
        T    = Tl(ind);
        p    = pl(ind);
        dT   = reshape(dTl(:, ind, :), [4, 1]);
        dp   = reshape(dpl(:, ind, :), [4, 1]); 
        beta = betal(ind);
    else
        rho = rhor(ind);
        v   = vr(ind);   
        T   = Tr(ind);
        p   = pr(ind);
        dT  = reshape(dTr(:, ind, :), [4, 1]);
        dp  = reshape(dpr(:, ind, :), [4, 1]);  
        beta = betar(ind);
    end
    
    % derivative of beta
    tmp      = 4 * (3 + tc.Cv * (2 * rho - 3)) ./ (tc.Cv * (rho - 3).^2);
    dbeta    = 2 * (dp ./ rho + tmp .* dT);
    dbeta(1) = dbeta(1) - 2 * (p ./ rho.^2 + 8 / tc.Cv * T .* (tc.Cv * rho + 3) ./ (rho - 3).^3);
    
    % derivative of global dissipation coefficient with respect to left and right state
    if (beta >= 0)
        sig = sign(v); sig(sig == 0) = 1;
        dalpha        = 0.5 ./ sqrt(beta) .* dbeta; 
        dalpha(1 : 2) = dalpha(1 : 2) + sig ./ rho .* [-v; 1];
    else
        tmp           = 0.5 ./ sqrt(v.^2 - beta);
        dalpha        = -tmp .* dbeta;
        dalpha(1 : 2) = dalpha(1 : 2) + tmp .* 2 .* v ./ rho .* [-v; 1];
    end    
    if (alphal(ind) >= alphar(ind))
        dalphal = dalpha;
        dalphar = zeros(4, 1);
    else     
        dalphal = zeros(4, 1);
        dalphar = dalpha;
    end   
elseif (inviscidFlux == 2)
    alphalIsChosen = alphal >= alphar;
    alpharIsChosen = ~alphalIsChosen;
    betalIsPositive = betal >= 0;
    betarIsPositive = betar >= 0;

    ind1 = alphalIsChosen & betalIsPositive;
    ind2 = alphalIsChosen & (~betalIsPositive);
    ind3 = alpharIsChosen & betarIsPositive;
    ind4 = alpharIsChosen & (~betarIsPositive);
    nf1  = sum(ind1);
    nf2  = sum(ind2);
    nf3  = sum(ind3);
    nf4  = sum(ind4);

    % derivative of beta
    tmpl = 4 * (3 + tc.Cv * (2 * rhol - 3)) ./ (tc.Cv * (rhol - 3).^2);
    tmpr = 4 * (3 + tc.Cv * (2 * rhor - 3)) ./ (tc.Cv * (rhor - 3).^2);
    dbetal = 2 * ([dpl(:, :, 1); dpl(:, :, 2); dpl(:, :, 3); dpl(:, :, 4)] ./ rhol + tmpl .* [dTl(:, :, 1); dTl(:, :, 2); dTl(:, :, 3); dTl(:, :, 4)]);
    dbetar = 2 * ([dpr(:, :, 1); dpr(:, :, 2); dpr(:, :, 3); dpr(:, :, 4)] ./ rhor + tmpr .* [dTr(:, :, 1); dTr(:, :, 2); dTr(:, :, 3); dTr(:, :, 4)]);
    dbetal(1, :) = dbetal(1, :) - 2 * (pl ./ rhol.^2 + 8 / tc.Cv * Tl .* (tc.Cv * rhol + 3) ./ (rhol - 3).^3);
    dbetar(1, :) = dbetar(1, :) - 2 * (pr ./ rhor.^2 + 8 / tc.Cv * Tr .* (tc.Cv * rhor + 3) ./ (rhor - 3).^3);

    % derivative of dissipation coefficient with respect to left and right state
    dalphal = zeros(4, nf);
    dalphar = zeros(4, nf);

    sig = sign(vl(ind1)); sig(sig == 0) = 1;
    dalphal(:, ind1)     = 0.5 ./ sqrt(betal(ind1)) .* dbetal(:, ind1); 
    dalphal(1 : 2, ind1) = dalphal(1 : 2, ind1) + sig ./ rhol(ind1) .* [-vl(ind1); ones(1, nf1)];

    tmpl = 0.5 ./ sqrt(vl(ind2).^2 - betal(ind2));
    dalphal(:, ind2)     = -tmpl .* dbetal(:, ind2);
    dalphal(1 : 2, ind2) = dalphal(1 : 2, ind2) + tmpl .* 2 .* vl(ind2) ./ rhol(ind2) .* [-vl(ind2); ones(1, nf2)];

    sig = sign(vr(ind3)); sig(sig == 0) = 1;
    dalphar(:, ind3)     = 0.5 ./ sqrt(betar(ind3)) .* dbetar(:, ind3); 
    dalphar(1 : 2, ind3) = dalphar(1 : 2, ind3) + sig ./ rhor(ind3) .* [-vr(ind3); ones(1, nf3)];

    tmpr = 0.5 ./ sqrt(vr(ind4).^2 - betar(ind4));
    dalphar(:, ind4)     = -tmpr .* dbetar(:, ind4);
    dalphar(1 : 2, ind4) = dalphar(1 : 2, ind4) + tmpr .* 2 .* vr(ind4) ./ rhor(ind4) .* [-vr(ind4); ones(1, nf4)];
end

end
