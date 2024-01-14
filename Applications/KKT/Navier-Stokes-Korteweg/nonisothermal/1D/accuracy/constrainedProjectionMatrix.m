function cpm = constrainedProjectionMatrix(msh, U, Lamb, dU, dLamb, rnorm, A, tc, bs, KKT, newton, mtol)

% row, column and value indices for the mass matrix
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);
ME = mychop(bs.phitw{1} * bs.phi{1}, 1.e-14);
r11 = re1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
c11 = ce1 + (0 : 3 * msh.nLElems - 1) * bs.nb;
v11 = repmat(msh.elemJac(:, msh.LElems) .* ME(:), [1, 3]);

% Add positivity constraints contribution
if (KKT.positivity)        
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    
    % Evaluate physical quantities
    u  = bs.phi{1} * reshape(U, [bs.nb, 3 * msh.nLElems]);
    w  = bs.phi{1} * reshape(A * U(1 : end / 3), [bs.nb, msh.nLElems]);
    v  = u(:, msh.nLElems + 1 : 2 * msh.nLElems) ./ u(:, 1 : msh.nLElems);
    T  = getTemperature(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, tc);
    p  = getPressure(u(:, 1 : msh.nLElems), T);
    dp = getPressureDeriv(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, tc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);  
    % pre-calculate some matrices
    rtmp = re1 + (0 : msh.nLElems - 1) * bs.nb;
    ctmp = ce1 + (0 : msh.nLElems - 1) * bs.nb; 
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, 1, 4, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp11 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, 4, 4, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp12 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    Mtmp1 = Mtmp11 * A + A' * (Mtmp11 + Mtmp12 * A);
    [rKKT1, cKKT1, vKKT1] = find(Mtmp1);
    
    rtmp = re2 + (0 : msh.nLElems - 1) * bs.nep;
    ctmp = ce2 + (0 : msh.nLElems - 1) * bs.nb;  
    vtmp = repelem(dp(:, :, 4), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]);
    
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
    
    % block-(1,1)
    for n = [1, 2, 3, 5, 6, 9]
        i = ceil(n / 4); j = n - (i - 1) * 4;
        dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v, u(:, 2 * msh.nLElems + 1 : end), w, i, j, tc);

        s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
        rKKT(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        cKKT(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
    end      
    
    % KKT-block (3,1)
    s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
    stmp = s;
    rKKT(s : e) = 3 * KKT.NDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
    cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
    vKKT(s : e) = repmat(bs.phi{1}(:), [msh.nLElems, 1]);
    for j = 1 : 3      
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 3 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp(:, :, j), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 4 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    [~, set2] = getActiveSets(u(:, 1 : msh.nLElems), p, Lamb, dU, dLamb, [rKKT(stmp : etmp); 3 * KKT.NDOF + KKT.LDOF + rKKT2], [cKKT(stmp : etmp); cKKT2], [vKKT(stmp : etmp); vKKT2], KKT, newton);    
      
    % KKT-block (3,3)
    s = e + 1; l = 2 * KKT.LDOF; e = e + l;
    rKKT(s : e) = 3 * KKT.NDOF + (1 : 2 * KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    vKKT = mychop(vKKT, mtol);
    cpm.Mat = sparse([r11(:); rKKT], [c11(:); cKKT], [v11(:); vKKT], 3 * KKT.NDOF + 2 * KKT.LDOF, 3 * KKT.NDOF + 2 * KKT.LDOF);
    cpm.Mat([false(3 * KKT.NDOF, 1); set2], 1 : 3 * KKT.NDOF) = 0;
    
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix 
        cpm.GT  = cpm.Mat';
        cpm.Mat = cpm.GT * cpm.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(cpm.Mat, 1));   
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        cpm.GT = speye(size(cpm.Mat, 1));
    end
else
    cpm.GT  = speye(3 * KKT.NDOF);
    cpm.Mat = sparse(r11, c11, v11, 3 * KKT.NDOF, 3 * KKT.NDOF);   
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

