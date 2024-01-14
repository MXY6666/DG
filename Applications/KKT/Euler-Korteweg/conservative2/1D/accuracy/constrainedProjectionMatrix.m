function cpm = constrainedProjectionMatrix(msh, U, Lamb, dU, dLamb, rnorm, bs, KKT, newton)

% row, column and value indices for the mass matrix
EM = computeElementMatrix_refLine(0, 0, bs, 1.e-16);
[re1, ce1, ve1] = find(EM);
r11 = re1 + (0 : 2 * msh.nLElems - 1) * bs.nb;
c11 = ce1 + (0 : 2 * msh.nLElems - 1) * bs.nb;
v11 = repmat(msh.elemJac(:, msh.LElems) .* ve1, [1, 2]);

% Add positivity constraints contribution
if (KKT.positivity)    
    % Get active sets 
    rho = bs.phi{1} * reshape(U(1 : KKT.NDOF), [bs.nb, msh.nLElems]);
    [~, set2] = getActiveSets(msh, rho, Lamb, dU, dLamb, 1, bs, KKT, newton);
    
    % Initialize row, column and value array for the KKT contributions
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);
    rKKT = zeros(2 * KKT.LDOF * bs.nb + KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;
    
    % KKT-block (3,1)
    s = 1; l = KKT.LDOF * bs.nb; e = l;
    stmp = s; etmp = e;
    rKKT(s : e) = 2 * KKT.NDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
    cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
    vKKT(s : e) = repmat(bs.phi{1}, [1, msh.nLElems]);    
    
    % KKT-block (1,3)
    s = e + 1; l = KKT.LDOF * bs.nb; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % KKT-block (3,3)
    s = e + 1; l = KKT.LDOF; e = e + l;
    rKKT(s : e) = 2 * KKT.NDOF + (1 : KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    cpm.Mat = sparse([r11(:); rKKT], [c11(:); cKKT], [v11(:); vKKT], 2 * KKT.NDOF + KKT.LDOF, 2 * KKT.NDOF + KKT.LDOF);
    cpm.Mat([false(2 * KKT.NDOF, 1); set2], 1 : 2 * KKT.NDOF) = 0;
    
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix 
        cpm.GT  = cpm.Mat';
        cpm.Mat = cpm.GT * cpm.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(cpm.Mat, 1));   
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        cpm.GT = speye(size(cpm.Mat, 1));
    end
else
    cpm.GT  = speye(2 * KKT.NDOF);
    cpm.Mat = sparse(r11, c11, v11, 2 * KKT.NDOF, 2 * KKT.NDOF);   
end

% Apply row-column scaling
if (newton.scaleMat)
    [cpm.Mat, cpm.PR, cpm.PC] = rowcolScale(cpm.Mat, 1.e-2, 5);
else
    si = size(cpm.Mat, 1);
    cpm.PR = speye(si);
    cpm.PC = speye(si);
end

end

