function Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, bs, bt, KKT, newton, mtol)

% Compute DG coefficients for current Runge-Kutta stage
U = getStageSolution(U0, KStage, dt, stage, bt);
U = reshape(U, [bs.nb, msh.nLElems]);

% Initialize the row, column and value array 
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);
r11 = zeros((2 * msh.nLElems + 4 * md.nIntLFaces) * bs.nb2, 1);
c11 = r11;
v11 = r11;

%**************************************************************************
%                       element contributions   
%**************************************************************************
u  = bs.phi{1} * U;

s = 1; l = msh.nLElems * bs.nb2; e = l;
r11(s : e) = re1 + (0 : msh.nLElems - 1) * bs.nb;
c11(s : e) = ce1 + (0 : msh.nLElems - 1) * bs.nb;
v11(s : e) = -bs.phitw{2} * (repelem(u, 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));

%**************************************************************************
%                        face contributions   
%**************************************************************************
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = bs.phi_face{1, 2} * U(:, leLIDs); 
ur = bs.phi_face{1, 1} * U(:, reLIDs);

% dissipation coefficient and its derivative 
[alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul, ur);

H_ul = 0.5 * (ul + alpha + dalphal .* (ul - ur));
H_ur = 0.5 * (ur - alpha + dalphar .* (ul - ur));

s = e + 1; l = 4 * nf * bs.nb2; e = e + l;
r11(s : e) = re1 + [leLIDs - 1, leLIDs - 1, reLIDs - 1, reLIDs - 1] * bs.nb;   
c11(s : e) = ce1 + [leLIDs - 1, reLIDs - 1, leLIDs - 1, reLIDs - 1] * bs.nb;
v11(s : e) = [bs.phitw_face{1, 2} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
              bs.phitw_face{1, 2} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), ...
             -bs.phitw_face{1, 1} * (repelem(H_ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), ...
             -bs.phitw_face{1, 1} * (repelem(H_ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]))];
                    
v11(1 : e) = dt * bt.diagCoe * v11(1 : e); 
%**************************************************************************
%                        mass matrix contributions   
%**************************************************************************
s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
r11(s : e) = re1 + (0 : msh.nLElems - 1) * bs.nb;
c11(s : e) = ce1 + (0 : msh.nLElems - 1) * bs.nb;
v11(s : e) = msh.elemJac(:, msh.LElems) .* reshape(bs.phitw{1} * bs.phi{1}, [bs.nb2, 1]);
v11 = mychop(v11, mtol);

%**************************************************************************
%                       KKT constraints contribution 
%**************************************************************************
% Add contribution of constraint equations
if (KKT.enable)
    nnep = (KKT.positivity + KKT.maximum) * bs.nep;
    if (bt.form == 1)
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;
    end
    
    % Initialize row, column and value array for the KKT contributions
    re2  = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2  = repelem((1 : bs.nb)', bs.nep);
    rKKT = zeros(2 * KKT.LDOF * bs.nb + KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;
        
    % KKT-block (3,1)
    e = 0; 
    if (KKT.positivity)
        s = e + 1; l = msh.nLElems * bs.nep * bs.nb; e = e + l;
        
        rKKT(s : e) = KKT.NDOF + re2 + (0 : msh.nLElems - 1) * nnep;
        cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repmat(cdt * bs.phi{1}, [1, msh.nLElems]); 
    end
    if (KKT.maximum)
        s = e + 1; l = msh.nLElems * bs.nep * bs.nb; e = e + l;
        
        rKKT(s : e) = KKT.NDOF + KKT.positivity * bs.nep + re2 + (0 : msh.nLElems - 1) * nnep;
        cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repmat(-cdt * bs.phi{1}, [1, msh.nLElems]); 
    end   
    stmp = 1; etmp = e;  
    
    % KKT-block (1,3)
    s = e + 1; l = KKT.LDOF * bs.nb; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets 
    [~, set2] = getActiveSets(msh, u, Lamb, dU, dLamb, cdt, bs, KKT, newton);
    
    % KKT-block (3,3)
    s = e + 1; l = KKT.LDOF; e = e + l;
    rKKT(s : e) = KKT.NDOF + (1 : KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    Imp.Mat = sparse([r11; rKKT], [c11; cKKT], [v11; vKKT], KKT.NDOF + KKT.LDOF, KKT.NDOF + KKT.LDOF);
    Imp.Mat([false(KKT.NDOF, 1); set2], 1 : KKT.NDOF) = 0;
     
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix
        Imp.GT  = Imp.Mat';
        Imp.Mat = Imp.GT * Imp.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(Imp.Mat, 1));
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        Imp.GT = speye(size(Imp.Mat, 1));
    end
else
    Imp.GT  = speye(KKT.NDOF);
    Imp.Mat = sparse(r11, c11, v11, KKT.NDOF, KKT.NDOF);
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
%                            subroutine   
%************************************************************************** 
% Compute dissipation coefficient and its derivatives
function [alpha, dalphal, dalphar] = getDissCoeAndDeriv(ul, ur)

% number of faces
nf = length(ul);

% dissipation coefficient
alphal = abs(ul);
alphar = abs(ur);
alpha  = max(alphal, alphar);

% derivative of dissipation coefficient with respect to left and right state
dalphal       = zeros(1, nf);
dalphar       = zeros(1, nf);
ind           = alphal >= alphar;
dalphal(ind)  = sign(ul(ind));
dalphar(~ind) = sign(ur(~ind));

end


