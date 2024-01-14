function [U, Lamb, fail] = constrainedProjection(msh, U, nt, t, flag, tc, quad2, bs, KKT, newton)

fprintf('\nConstrained projection\n\n')
% Initialize Newton solution vector
if (KKT.enable)
    Lamb = zeros(KKT.LDOF, 1);
else
    Lamb = [];
end
V = zeros(KKT.NDOF + KKT.LDOF, 1);
V(1 : KKT.NDOF) = U; 

% Newton step
fail  = true;
U0    = U;
res   = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad2, bs, KKT);
norm0 = norm(res);
rnorm = 1;
dV    = zeros(KKT.NDOF + KKT.LDOF, 1);
for k = 1 : newton.maxit
    if (k == 1 || L2_res > min(newton.tolres, 1.e-12))
        % Compute the Jacobian-associated matrix 
        cpm = constrainedProjectionMatrix(msh, U, Lamb, dV(1 : KKT.NDOF), dV(KKT.NDOF + 1 : end), rnorm, bs, KKT, newton);

        % Multiply for least-squares method the residual vector with the transposed of the original implicit matrix
        resLS = cpm.GT * res;

        % Scale residual vectors
        scaRes   = cpm.PR * res;
        scaResLS = cpm.PR * resLS;        
        
        % Set switch to indicate that full Newton is used
        fullNewton = true;         
    else
        % Multiply for least-squares method the residual vector with the transposed of the original implicit matrix
        resLS = cpm.GT * res;
        
        % Scale residual vectors
        scaResLS = cpm.PR * resLS;  
        
        % Set switch to indicate that full Newton is used
        fullNewton = false;            
    end
    
    if (k == 1)
        % Compute merit function
        theta = 0.5 * (scaRes' * scaRes);
        
        % Initialize merit function reference values
        theta_refs = theta * ones(newton.nref, 1);
        
        % Print convergence history 
        fprintf('Time-step      Time       Newton-iter   FN    L2-norm Residual   L2-norm DelV\n')  
        fprintf('%6d      %10.4e      %3d       %4d    %13.4e\n', nt, t, 0, fullNewton, snorm(scaRes))         
    end
    
    % Solve linear system for Newton correction 
    scadV = linSolve(cpm.Mat, -scaResLS, newton.linearSolver); 
    dV    = cpm.PC * scadV;
    
    % Store solution and residual vector of previous Newton iteration
    Vold      = V;
    scaResold = scaRes;
    
    % Compute trial solution
    V = Vold + dV;
    U = V(1 : KKT.NDOF);
    Lamb = max(V(KKT.NDOF + 1 : end), 0);
    V(KKT.NDOF + 1 : end) = Lamb;
      
    % Compute residual vector for new descent vector
    res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad2, bs, KKT);
    scaRes = cpm.PR * res;
    
    % Compute merit function for new descent vector
    theta = 0.5 * (scaRes' * scaRes);
        
    testnorm = norm(scaResold, Inf);
    if(norm(scaRes, Inf) >= newton.gamma * norm(scaResold, Inf) && testnorm > 1.e-14)  
        % Get reference merit function
        theta_ref = max(theta_refs);
                
        if (theta - theta_ref > -newton.sigma * (scadV' * scadV) && theta >= 1.e-20)
            for l = 1 : newton.maxitLineSearch - 1
                % Compute trial solution
                V = Vold + newton.beta^l * dV;
                U = V(1 : KKT.NDOF);
                Lamb = max(V(KKT.NDOF + 1 : end), 0);
                V(KKT.NDOF + 1 : end) = Lamb;

                % Compute residual vector for new descent vector
                res = constrainedProjectionResidual(msh, U0, U, Lamb, flag, tc, quad2, bs, KKT);
                scaRes = cpm.PR * res;        

                % Compute merit function for new descent vector
                theta = 0.5 * (scaRes' * scaRes);   

                if (theta - theta_ref <= -newton.sigma * newton.beta^l * (scadV' * scadV) || theta < 1.e-20)
                    % Found proper search direction
                    break               
                end

                if (l == newton.maxitLineSearch - 1)
                    % Line search did not converge
                    fprintf('\nNo feasible search direction found in constrained projection\n')
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
    
    % relative norm of function vector
    rnorm = norm(res) / (norm0 + eps); 
    
    % Compute L2-norm of residual and solution difference
    L2_res = snorm(scaRes);
    L2_dV  = snorm(scadV);  
    
    % Print convergence history 
    fprintf('%6d      %10.4e      %3d       %4d    %13.4e      %11.4e\n', nt, t, k, fullNewton, L2_res, L2_dV)
    
    % Determine if stop the Newton method
    if (L2_res <= newton.tolres && L2_dV <= newton.tolres)
        fail = false;
        return
    end    
end
fprintf('\nNewton method is converging, but did not reduce the error to the requested level\n')
fprintf('\nIncrease number of Newton interations\n')

end