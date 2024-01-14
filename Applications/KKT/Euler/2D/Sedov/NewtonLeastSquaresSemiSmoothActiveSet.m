function [KStage, Lamb, LStage, k, fail] = NewtonLeastSquaresSemiSmoothActiveSet(msh, md, U0, KStage, Lamb, LStage, nt, t, dt, stage, M, tc, bs, bt, KKT, newton, mtol)

% Ensure that initial solution for active set Newton method is realizable
Lamb = checkCompatibility(msh, U0, KStage, Lamb, dt, stage, tc, bs, bt, KKT);

% Initialize Newton solution vector
V = zeros(4 * KKT.NDOF + 2 * KKT.LDOF, 1);
V(1 : 4 * KKT.NDOF) = KStage(:, stage.Ks); 
V(4 * KKT.NDOF + 1 : end) = Lamb;

% Newton step
[res, ~, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, tc, bs, bt, KKT);
norm0 = norm(res);
rnorm = 1;
dV    = zeros(4 * KKT.NDOF + 2 * KKT.LDOF, 1);
fail  = true;
for k = 1 : newton.maxit
    if (k == 1 || L2_res > min(newton.tolres, 1.e-12))
        % Compute the Jacobian-associated matrix 
        Imp = computeImplicitMatrix(msh, md, U0, KStage, Lamb, dV(1 : 4 * KKT.NDOF), dV(4 * KKT.NDOF + 1 : end), dt, stage, rnorm, tc, bs, bt, KKT, newton, mtol);

        % Multiply for least-squares method the residual vector with the transposed of the original implicit matrix
        resLS = Imp.GT * res;

        % Scale residual vectors
        scaRes   = Imp.PR * res;
        scaResLS = Imp.PR * resLS;   
        
        % Set switch to indicate that full Newton is used
        fullNewton = true;
    else
        % Multiply for least-squares method the residual vector with the transposed of the original implicit matrix
        resLS = Imp.GT * res;
        
        % Scale residual vectors
        scaResLS = Imp.PR * resLS;
        
        % Set switch to indicate that full Newton is used
        fullNewton = false;
    end
    
    if (k == 1)
        % Compute merit function
        theta = 0.5 * (scaRes' * scaRes);
        
        % Initialize merit function reference values
        theta_refs = theta * ones(newton.nref, 1);
        
        % Compute CFL number
        CFL = waveVel * dt / md.meshSize;
        
        % Print convergence history 
        fprintf('\n\nTime-step      Time       RK-stage     Newton-iter      CFL       FN    L2-norm Residual   L2-norm DelV\n')  
        fprintf('%6d      %10.4e      %3d           %3d       %10.4e   %3d      %12.4e\n', nt, t + dt, stage.s, 0, CFL, fullNewton, snorm(scaRes))
    end
    
    % Solve linear system for Newton correction 
    scadV = linSolve(Imp.Mat, -scaResLS, newton.linearSolver); 
    dV    = Imp.PC * scadV;
    
    % Store solution and residual vector of previous Newton iteration
    Vold      = V;
    scaResold = scaRes;
    
    % Compute trial solution
    V = Vold + dV;
    KStage(:, stage.Ks) = V(1 : 4 * KKT.NDOF);
    Lamb = max(V(4 * KKT.NDOF + 1 : end), 0);
    V(4 * KKT.NDOF + 1 : end) = Lamb;
      
    % Compute residual vector for new descent vector
    [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, tc, bs, bt, KKT);
    scaRes = Imp.PR * res;

    % Compute merit function for new descent vector
    theta = 0.5 * (scaRes' * scaRes);  
    
    if (norm(scaRes, Inf) >= newton.gamma * norm(scaResold, Inf))  
        % Get reference merit function
        theta_ref = max(theta_refs);   
        
        % Perform line search
        if (theta - theta_ref > -newton.sigma * (scadV' * scadV) && theta >= 1.e-15)
            for l = 1 : newton.maxitLineSearch - 1
                % Compute trial solution
                V = Vold + newton.beta^l * dV;
                KStage(:, stage.Ks) = V(1 : 4 * KKT.NDOF);
                Lamb = max(V(4 * KKT.NDOF + 1 : end), 0);
                V(4 * KKT.NDOF + 1 : end) = Lamb;

                % Compute residual vector for new descent vector
                [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, tc, bs, bt, KKT);
                scaRes = Imp.PR * res;        

                % Compute merit function for new descent vector
                theta = 0.5 * (scaRes' * scaRes);   

                if (theta - theta_ref <= -newton.sigma * newton.beta^l * (scadV' * scadV) || theta < 1.e-15)
                    % Found proper search direction
                    break               
                end

                if (l == newton.maxitLineSearch - 1)
                    % Line search did not converge
                    fprintf('\nNo feasible search direction found  -- restart Newton iteration with small time step\n')
                    return
                end
            end   
        end
    end
       
    % Update reference merit function vector
    [theta_ref, ind] = max(theta_refs);
    if (theta < theta_ref)
        theta_refs(ind) = theta;
    end
    
    % relative norm of residual vector
    rnorm = norm(res) / (norm0 + eps);    
    
    % Compute CFL number
    CFL = waveVel * dt / md.meshSize;
    
    % Compute L2-norm of residual and solution difference
    L2_res = snorm(scaRes);
    L2_dV  = snorm(scadV);
    
    % Print convergence history 
    fprintf('%6d      %10.4e      %3d           %3d       %10.4e   %3d      %12.4e    %12.4e\n', nt, t + dt, stage.s, k, CFL, fullNewton, L2_res, L2_dV)   
      
    % Determine if stop the Newton method
    if (L2_res <= newton.tolres && L2_dV <= newton.toldV)
        fail = false;
        return
    end
    
    % Restart Newton iteration if a local minimum is found
    if (L2_res > newton.tolres && L2_dV <= newton.toldVLocal)
        fprintf('\nLocal minimum found  -- restart Newton iteration with small time step\n')  
        return
    end    
end
fprintf('\nNewton method is converging, but did not reduce the error to the requested level\n')
fprintf('\nIncrease number of Newton interations\n')

end