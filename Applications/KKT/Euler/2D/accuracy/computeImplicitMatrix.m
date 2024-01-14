function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, tc, bs, bt, KKT, newton, mtol)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);
r11 = zeros((16 * (msh.nLElems + 4 * md.nIntLFaces) + 4 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

%**************************************************************************
%                       element contributions   
%**************************************************************************
u = bs.phi{1} * U;
v = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];

e = 0;
for n = 1 : 16
    if (n == 1 || n == 4)
        continue;
    end
    i = ceil(n / 4); j = n - (i - 1) * 4;
    ent1 = getJacobianEntry1(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), i, j, tc);
    ent2 = getJacobianEntry2(u(:, 1 : 4 : end), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 4 : 4 : end), i, j, tc);
    
    s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
    r11(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s : e) = -bs.phitw{2} * (repelem(msh.elemJxix( :, msh.LElems) .* ent1, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])) + ...
                 -bs.phitw{3} * (repelem(msh.elemJetay(:, msh.LElems) .* ent2, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
end

%**************************************************************************
%                        face contributions   
%**************************************************************************
% internal vertical face contributions
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
vr = ur(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul(1 : bs.nfp, :), ur(1 : bs.nfp, :), vl(1 : bs.nfp, :), vr(1 : bs.nfp, :), vl(bs.nfp + 1 : 2 * bs.nfp, :), vr(bs.nfp + 1 : 2 * bs.nfp, :), pl, pr, 1, tc);
                 
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    
    entl = getJacobianEntry1(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), i, j, tc);
    entr = getJacobianEntry1(ur(1 : bs.nfp, :), vr(1 : bs.nfp, :), vr(bs.nfp + 1 : end, :), ur(3 * bs.nfp + 1 : end, :), i, j, tc);   
    H_ul = 0.5 * (entl - dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s = e + 1; l = 4 * nf * bs.nb2; e = e + l;
    r11(s : e) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    c11(s : e) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb; 
    v11(s : e) = [bs.phitw_face{1, 2} * (repelem(J .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                  bs.phitw_face{1, 2} * (repelem(J .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
                 -bs.phitw_face{1, 1} * (repelem(J .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                 -bs.phitw_face{1, 1} * (repelem(J .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];    
end

% internal horizontal face contributions
faceIDs = md.intLFaces{2, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs), [4 * bs.nfp, nf]);
vl = ul(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ul(1 : bs.nfp, :), [2, 1]);
vr = ur(bs.nfp + 1 : 3 * bs.nfp, :) ./ repmat(ur(1 : bs.nfp, :), [2, 1]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul(1 : bs.nfp, :), ur(1 : bs.nfp, :), vl(1 : bs.nfp, :), vr(1 : bs.nfp, :), vl(bs.nfp + 1 : 2 * bs.nfp, :), vr(bs.nfp + 1 : 2 * bs.nfp, :), pl, pr, 2, tc);
                 
for n = 1 : 16
    i = ceil(n / 4); j = n - (i - 1) * 4;
    
    entl = getJacobianEntry2(ul(1 : bs.nfp, :), vl(1 : bs.nfp, :), vl(bs.nfp + 1 : end, :), ul(3 * bs.nfp + 1 : end, :), i, j, tc);
    entr = getJacobianEntry2(ur(1 : bs.nfp, :), vr(1 : bs.nfp, :), vr(bs.nfp + 1 : end, :), ur(3 * bs.nfp + 1 : end, :), i, j, tc);   
    H_ul = 0.5 * (entl - dalphal((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar((j - 1) * bs.nfp + 1 : j * bs.nfp, :) .* (ur((i - 1) * bs.nfp + 1 : i * bs.nfp, :) - ul((i - 1) * bs.nfp + 1 : i * bs.nfp, :)) - alpha * (i == j));
    
    s = e + 1; l = 4 * nf * bs.nb2; e = e + l;
    r11(s : e) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    c11(s : e) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb; 
    v11(s : e) = [bs.phitw_face{1, 4} * (repelem(J .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf])), ...
                  bs.phitw_face{1, 4} * (repelem(J .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf])), ...
                 -bs.phitw_face{1, 3} * (repelem(J .* H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 4}, [1, nf])), ...
                 -bs.phitw_face{1, 3} * (repelem(J .* H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 3}, [1, nf]))];    
end
                    
v11(1 : e) = dt * bt.diagCoe * v11(1 : e); 
%**************************************************************************
%                        mass matrix contributions   
%**************************************************************************
s = e + 1; l = 4 * msh.nLElems * bs.nb2; e = e + l;
r11(s : e) = re1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
c11(s : e) = ce1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
v11(s : e) = repmat(msh.elemJac(:, msh.LElems) .* reshape(bs.phitw{1} * bs.phi{1}, [bs.nb2, 1]), [1, 4]);

r11(e + 1 : end) = [];
c11(e + 1 : end) = [];
v11(e + 1 : end) = [];
v11 = mychop(v11, mtol);

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
    vKKT = mychop(vKKT, mtol);
    Imp.Mat = sparse([r11; rKKT], [c11; cKKT], [v11; vKKT], 4 * KKT.NDOF + 2 * KKT.LDOF, 4 * KKT.NDOF + 2 * KKT.LDOF);
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
    Imp.Mat = sparse(r11, c11, v11, 4 * KKT.NDOF, 4 * KKT.NDOF);    
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
%                                subroutine   
%************************************************************************** 
% Compute dissipation coefficient and its derivatives at internal faces
function [alpha, dalphal, dalphar] = getDissCoeAndDeriv(rhol, rhor, v1l, v1r, v2l, v2r, pl, pr, dir, tc)

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

dalphal = repmat(0.5 * sqrt(tc.gamma ./ rhol ./ pl), [4, 1]) .* [0.5 * (tc.gamma - 1) * (v1l.^2 + v2l.^2) - pl ./ rhol; (1 - tc.gamma) * v1l; (1 - tc.gamma) * v2l; repmat(tc.gamma - 1, [nfp, nf_])];
dalphar = repmat(0.5 * sqrt(tc.gamma ./ rhor ./ pr), [4, 1]) .* [0.5 * (tc.gamma - 1) * (v1r.^2 + v2r.^2) - pr ./ rhor; (1 - tc.gamma) * v1r; (1 - tc.gamma) * v2r; repmat(tc.gamma - 1, [nfp, nf_])];

tmp = sigl ./ rhol;
dalphal(1 : nfp, :) = dalphal(1 : nfp, :) - tmp .* vnl;
dalphal(ind, :)     = dalphal(ind, :)     + tmp;

tmp = sigr ./ rhor;
dalphar(1 : nfp, :) = dalphar(1 : nfp, :) - tmp .* vnr;
dalphar(ind, :)     = dalphar(ind, :)     + tmp;

dalphal(repmat(indr, [4, 1])) = 0;
dalphar(repmat(indl, [4, 1])) = 0;

end
