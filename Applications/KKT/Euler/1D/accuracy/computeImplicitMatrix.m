function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, tc, bs, bt, KKT, newton, mtol)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);

r11 = zeros(((msh.nLElems + 4 * md.nIntLFaces) * 9 + 3 * msh.nLElems) * bs.nb2, 1);
c11 = r11;
v11 = r11;

%**************************************************************************
%                       element contributions   
%**************************************************************************
u = bs.phi{1} * U;
v = u(:, 2 : 3 : end) ./ u(:, 1 : 3 : end);

e = 0;
for n = 1 : 9
    if (n == 1 || n == 3)
        continue;
    end
    i = ceil(n / 3); j = n - (i - 1) * 3;
    ent = getJacobianEntry(u(:, 1 : 3 : end), v, u(:, 3 : 3 : end), i, j, tc);
    
    s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
    r11(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
    c11(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
    v11(s : e) = -bs.phitw{2} * (repelem(ent, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
end

%**************************************************************************
%                        face contributions   
%**************************************************************************
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

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul(1, :), ur(1, :), vl, vr, pl, pr, tc);
                 
for n = 1 : 9
    i = ceil(n / 3); j = n - (i - 1) * 3;
    
    entl = getJacobianEntry(ul(1, :), vl, ul(3, :), i, j, tc);
    entr = getJacobianEntry(ur(1, :), vr, ur(3, :), i, j, tc);
    H_ul = 0.5 * (entl - dalphal(j, :) .* (ur(i, :) - ul(i, :)) + alpha * (i == j));
    H_ur = 0.5 * (entr - dalphar(j, :) .* (ur(i, :) - ul(i, :)) - alpha * (i == j));
    
    s = e + 1; l = 4 * nf * bs.nb2; e = e + l;
    r11(s : e) = re1 + [(i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + leLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1, (i - 1) * msh.nLElems + reLIDs - 1] * bs.nb;
    c11(s : e) = ce1 + [(j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1, (j - 1) * msh.nLElems + leLIDs - 1, (j - 1) * msh.nLElems + reLIDs - 1] * bs.nb; 
    v11(s : e) = [bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                  bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
                 -bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
                 -bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];    
end
                    
v11(1 : e) = dt * bt.diagCoe * v11(1 : e); 
%**************************************************************************
%                        mass matrix contributions   
%**************************************************************************
s = e + 1; l = 3 * msh.nLElems * bs.nb2; e = e + l;
r11(s : e) = re1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
c11(s : e) = ce1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
v11(s : e) = repmat(msh.elemJac(:, msh.LElems) .* reshape(bs.phitw{1} * bs.phi{1}, [bs.nb2, 1]), [1, 3]);

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
    p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
        
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
    [~, set2] = getActiveSets(u(:, 1 : 3 : end), p, Lamb, dU, dLamb, rKKT(stmp : etmp), cKKT(stmp : etmp), vKKT(stmp : etmp), KKT, newton);    
      
    % KKT-block (3,3)
    s = e + 1; l = 2 * KKT.LDOF; e = e + l;
    rKKT(s : e) = 3 * KKT.NDOF + (1 : 2 * KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    vKKT = mychop(vKKT, mtol);
    Imp.Mat = sparse([r11; rKKT], [c11; cKKT], [v11; vKKT], 3 * KKT.NDOF + 2 * KKT.LDOF, 3 * KKT.NDOF + 2 * KKT.LDOF);
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
    Imp.Mat = sparse(r11, c11, v11, 3 * KKT.NDOF, 3 * KKT.NDOF);    
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
function [alpha, dalphal, dalphar] = getDissCoeAndDeriv(rhol, rhor, vl, vr, pl, pr, tc)
    
alphal = getEigenMax(rhol, vl, pl, tc);
alphar = getEigenMax(rhor, vr, pr, tc);
alpha  = max(alphal, alphar);

ind1 = alphal >= alphar;
ind2 = ~ind1;
nf1  = sum(ind1);
nf2  = sum(ind2);

% derivative of dissipation coefficient with respect to left and right state
nf = size(rhol, 2);
dalphal = zeros(3, nf);
dalphar = zeros(3, nf);

dalphal(:, ind1) = 0.5 * sqrt(tc.gamma ./ rhol(ind1) ./ pl(ind1)) .* [0.5 * (tc.gamma - 1) * vl(ind1).^2 - pl(ind1) ./ rhol(ind1); (1 - tc.gamma) * vl(ind1); repmat(tc.gamma - 1, [1, nf1])];
dalphar(:, ind2) = 0.5 * sqrt(tc.gamma ./ rhor(ind2) ./ pr(ind2)) .* [0.5 * (tc.gamma - 1) * vr(ind2).^2 - pr(ind2) ./ rhor(ind2); (1 - tc.gamma) * vr(ind2); repmat(tc.gamma - 1, [1, nf2])];

sig = sign(vl(ind1)); sig(sig == 0) = 1;
dalphal(1 : 2, ind1) = dalphal(1 : 2, ind1) + sig ./ rhol(ind1) .* [-vl(ind1); ones(1, nf1)];

sig = sign(vr(ind2)); sig(sig == 0) = 1;
dalphar(1 : 2, ind2) = dalphar(1 : 2, ind2) + sig ./ rhor(ind2) .* [-vr(ind2); ones(1, nf2)];

end

