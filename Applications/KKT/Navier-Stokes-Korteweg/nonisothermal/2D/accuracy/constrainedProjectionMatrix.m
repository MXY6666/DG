function cpm = constrainedProjectionMatrix(msh, U, Lamb, dU, dLamb, rnorm, Ax, Ay, tc, bs, KKT, newton, mtol)

% row, column and value indices for the mass matrix
re1 = repmat((1 : bs.nb)', [bs.nb, 1]);
ce1 = repelem((1 : bs.nb)', bs.nb);
ME = mychop(bs.phitw{1} * bs.phi{1}, 1.e-14);
r11 = re1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
c11 = ce1 + (0 : 4 * msh.nLElems - 1) * bs.nb;
v11 = repmat(msh.elemJac(:, msh.LElems) .* ME(:), [1, 4]);

% Add positivity constraints contribution
if (KKT.positivity)   
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    
    % Evaluate physical quantities
    u  = bs.phi{1} * reshape(U, [bs.nb, 4 * msh.nLElems]);
    w1 = bs.phi{1} * reshape(Ax * U(1 : end / 4), [bs.nb, msh.nLElems]);
    w2 = bs.phi{1} * reshape(Ay * U(1 : end / 4), [bs.nb, msh.nLElems]);    
    v  = u(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(u(:, 1 : msh.nLElems), [1, 2]);
    T  = getTemperature(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, tc);
    p  = getPressure(u(:, 1 : msh.nLElems), T);
    dp = getPressureDeriv(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, tc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    re2 = repmat((1 : bs.nep)', [bs.nb, 1]);
    ce2 = repelem((1 : bs.nb)', bs.nep);    
    
    % pre-calculate some matrices
    rtmp = re1 + (0 : msh.nLElems - 1) * bs.nb;
    ctmp = ce1 + (0 : msh.nLElems - 1) * bs.nb;
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, 1, 5, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp11 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, 1, 6, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp12 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, 5, 5, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp13 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, 6, 6, tc);
    vtmp = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));
    Mtmp14 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.NDOF, KKT.NDOF);
    
    Mtmp1 = Mtmp11 * Ax + Mtmp12 * Ay + Ax' * (Mtmp11 + Mtmp13 * Ax) + Ay' * (Mtmp12 + Mtmp14 * Ay);
    [rKKT1, cKKT1, vKKT1] = find(Mtmp1);
    
    rtmp  = re2 + (0 : msh.nLElems - 1) * bs.nep;
    ctmp  = ce2 + (0 : msh.nLElems - 1) * bs.nb;  
    
    vtmp = repelem(dp(:, :, 5), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]);
    Mtmp21 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.LDOF, KKT.NDOF);
    
    vtmp = repelem(dp(:, :, 6), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]);
    Mtmp22 = sparse(rtmp(:), ctmp(:), vtmp(:), KKT.LDOF, KKT.NDOF);
    
    Mtmp2 = Mtmp21 * Ax + Mtmp22 * Ay;
    [rKKT2, cKKT2, vKKT2] = find(Mtmp2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Initialize row, column and value array for the KKT contributions 
    rKKT = zeros(length(vKKT1) + 2 * length(vKKT2) + 9 * msh.nLElems * bs.nb2 + 10 * msh.nLElems * bs.nb * bs.nep + 2 * KKT.LDOF, 1);
    cKKT = rKKT;
    vKKT = rKKT;    
    
    e = length(vKKT1) + 2 * length(vKKT2);
    rKKT(1 : e) = [rKKT1; 4 * KKT.NDOF + KKT.LDOF + rKKT2; cKKT2];
    cKKT(1 : e) = [cKKT1; cKKT2; 4 * KKT.NDOF + KKT.LDOF + rKKT2];
    vKKT(1 : e) = [vKKT1; vKKT2; -vKKT2];  
       
    % block-(1,1)
    for n = [1, 2, 3, 4, 7, 8, 13, 15, 19]
        i = ceil(n / 6); j = n - (i - 1) * 6;
        dp2 = getPressureHessian(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, i, j, tc);

        s = e + 1; l = msh.nLElems * bs.nb2; e = e + l;
        rKKT(s : e) = re1 + ((i - 1) * msh.nLElems : i * msh.nLElems - 1) * bs.nb;
        cKKT(s : e) = ce1 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = -bs.phi{1}' * (repelem(dp2 .* Lamb(:, msh.nLElems + 1 : end), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]));  
    end      
    
    % KKT-block (3,1)
    s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
    stmp = s;
    rKKT(s : e) = 4 * KKT.NDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
    cKKT(s : e) = ce2 + (0 : msh.nLElems - 1) * bs.nb;
    vKKT(s : e) = repmat(bs.phi{1}(:), [msh.nLElems, 1]);
    for j = 1 : 4     
        s = e + 1; l = msh.nLElems * bs.nb * bs.nep; e = e + l;
        rKKT(s : e) = 4 * KKT.NDOF + KKT.LDOF + re2 + (0 : msh.nLElems - 1) * bs.nep;
        cKKT(s : e) = ce2 + ((j - 1) * msh.nLElems : j * msh.nLElems - 1) * bs.nb;
        vKKT(s : e) = repelem(dp(:, :, j), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems]);
    end
    etmp = e;
    
    % KKT-block (1,3)
    s = e + 1; l = 5 * msh.nLElems * bs.nb * bs.nep; e = e + l;
    rKKT(s : e) = cKKT(stmp : etmp);
    cKKT(s : e) = rKKT(stmp : etmp);
    vKKT(s : e) = -vKKT(stmp : etmp);
    
    % Get active sets
    [~, set2] = getActiveSets(u(:, 1 : msh.nLElems), p, Lamb, dU, dLamb, [rKKT(stmp : etmp); 4 * KKT.NDOF + KKT.LDOF + rKKT2], [cKKT(stmp : etmp); cKKT2], [vKKT(stmp : etmp); vKKT2], KKT, newton);    
      
    % KKT-block (3,3)
    s = e + 1; l = 2 * KKT.LDOF; e = e + l;
    rKKT(s : e) = 4 * KKT.NDOF + (1 : 2 * KKT.LDOF)';
    cKKT(s : e) = rKKT(s : e);
    vKKT(s : e) = set2(:);
    
    % Generate the Jacobian matrix
    vKKT = mychop(vKKT, mtol);
    cpm.Mat = sparse([r11(:); rKKT], [c11(:); cKKT], [v11(:); vKKT], 4 * KKT.NDOF + 2 * KKT.LDOF, 4 * KKT.NDOF + 2 * KKT.LDOF);
    cpm.Mat([false(4 * KKT.NDOF, 1); set2], 1 : 4 * KKT.NDOF) = 0;
    
    if (newton.linearSolver.type == 1 || newton.linearSolver.type == 3)
        % Compute least-squares normal matrix 
        cpm.GT  = cpm.Mat';
        cpm.Mat = cpm.GT * cpm.Mat + max(newton.alpha_bar * rnorm, 1.e-14) * speye(size(cpm.Mat, 1));   
    elseif (newton.linearSolver.type == 2)
        % Use lsqminnorm solver
        cpm.GT = speye(size(cpm.Mat, 1));
    end
else
    cpm.GT  = speye(4 * KKT.NDOF);
    cpm.Mat = sparse(r11, c11, v11, 4 * KKT.NDOF, 4 * KKT.NDOF);   
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

