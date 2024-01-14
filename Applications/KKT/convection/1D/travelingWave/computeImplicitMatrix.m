function Imp = computeImplicitMatrix(msh, U0, KStage, Lamb, dU, dLamb, dt, stage, rnorm, M, S, bs, bt, KKT, newton)

% DG discretization contribution  
M11 = M + dt * bt.diagCoe * S;

% Add contribution of constraint equations
if (KKT.enable)
    nnep = (KKT.positivity + KKT.maximum) * bs.nep;
    if (bt.form == 1)
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;
    end
    U = getStageSolution(U0, KStage, dt, stage, bt);
    u = bs.phi{1} * reshape(U, [bs.nb, msh.nLElems]);
        
    % Initialize row, column and value array for the KKT contributions
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);
    rKKT = zeros(2 * KKT.LDOF * bs.nb + KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;
    
    % block-(1,1)
    [r11, c11, v11] = find(M11);
    
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
    Imp.GT  = speye(size(M11, 1));
    Imp.Mat = M11;
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


