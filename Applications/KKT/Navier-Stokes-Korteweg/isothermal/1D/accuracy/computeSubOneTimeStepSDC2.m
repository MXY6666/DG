% Compute one time step in a time subinterval 
function [U, Lamb, AU, iL, L, it, fail] = computeSubOneTimeStepSDC2(msh, md, U0, Lamb, AU0, iL0, L0, k, m, nt, ts, dts, flag, flux, AE, A2, A3, A32, ME, IME, tc, quad1, bs, SDC, KKT, newton)
    
t = ts(1) + sum(dts);
% Ensure that initial solution for active set Newton method is realizable
Lamb = checkCompatibility(msh, U0, Lamb, bs, KKT);

% Initialize Newton solution vector
V = zeros(2 * KKT.NDOF + KKT.LDOF, 1);
V(1 : 2 * KKT.NDOF) = U0; 
V(2 * KKT.NDOF + 1 : end) = Lamb;

% Compute part of Jacobian matrix and residuals
if (k == 1)
    Jac1 = speye(2 * KKT.NDOF) - dts(m - 1) * AU0;
    res1 = -U0 - dts(m - 1) * iL0;    
else
    Jac1 = speye(2 * KKT.NDOF) - SDC.theta * dts(m - 1) * AU0;
    res1 = -U0 - SDC.theta * dts(m - 1) * (iL0 - L0(:, m)) - sum(dts) * sum(SDC.weights(m - 1, :) .* L0, 2); 
end

% Newton step
fail  = true;
res   = computeResidual(msh, md, U0, Lamb, AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
norm0 = norm(res);
rnorm = 1;
dV    = zeros(2 * KKT.NDOF + KKT.LDOF, 1);
for it = 1 : newton.maxit
    if (it == 1 || L2_res > min(newton.tolres, 1.e-12))
        % Compute the Jacobian-associated matrix 
        Imp = computeImplicitMatrix(msh, md, V(1 : 2 * KKT.NDOF), V(2 * KKT.NDOF + 1 : end), dV(1 : 2 * KKT.NDOF), dV(2 * KKT.NDOF + 1 : end), Jac1, k, m, dts, rnorm, flux, A3, A32, ME, IME, bs, SDC, KKT, newton);
        
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
    
    if (it == 1)        
        % Compute merit function
        theta = 0.5 * (scaRes' * scaRes);
        
        % Initialize merit function reference values
        theta_refs = theta * ones(newton.nref, 1);
        
        % Print convergence history 
        fprintf('\n\nTime-step      Time        k-m     Newton-iter   FN    L2-norm Residual   L2-norm DelU\n')  
        fprintf('%6d      %10.4e     %d-%d        %3d       %3d     %12.4e\n', nt, t, k, m - 1, 0, fullNewton, snorm(scaRes))
    end    
    
    % Solve linear system for Newton correction     
    scadV = linSolve(Imp.Mat, -scaResLS, newton.linearSolver); 
    dV    = Imp.PC * scadV;
    
    % Store solution and residual vector of previous Newton iteration
    Vold      = V;
    scaResold = scaRes;
    
    % Compute trial solution
    V = Vold + dV;
    V(2 * KKT.NDOF + 1 : end) = max(V(2 * KKT.NDOF + 1 : end), 0);
    
    % Compute residual vector for new descent vector
    [res, iL] = computeResidual(msh, md, V(1 : 2 * KKT.NDOF), V(2 * KKT.NDOF + 1 : end), AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
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
                V(2 * KKT.NDOF + 1 : end) = max(V(2 * KKT.NDOF + 1 : end), 0);

                % Compute residual vector for new descent vector
                [res, iL] = computeResidual(msh, md, V(1 : 2 * KKT.NDOF), V(2 * KKT.NDOF + 1 : end), AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
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
                    U    = V(1 : 2 * KKT.NDOF);
                    Lamb = V(2 * KKT.NDOF + 1 : end); 
                    AU   = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
                    iL   = zeros(2 * KKT.NDOF, 1); 
                    L    = zeros(2 * KKT.NDOF, 1);  
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
    
    % Compute L2-norm of residual and solution difference
    L2_res = snorm(scaRes);
    L2_dV  = snorm(scadV);
    
    % Print convergence history 
    fprintf('%6d      %10.4e     %d-%d        %3d       %3d     %12.4e      %12.4e\n', nt, t, k, m - 1, it, fullNewton, L2_res, L2_dV)
    
    % Determine if stop the Newton method
    if (L2_res <= newton.tolres && L2_dV <= newton.toldV)
        fail = false;
        U    = V(1 : 2 * KKT.NDOF);
        Lamb = V(2 * KKT.NDOF + 1 : end);
        L    = iL;      
        break
    end
    
    % Restart Newton iteration if a local minimum is found
    if (L2_res > newton.tolres && L2_dV <= newton.toldVLocal)
        fprintf('\nLocal minimum found  -- restart Newton iteration with small time step\n') 
        U    = V(1 : 2 * KKT.NDOF);
        Lamb = V(2 * KKT.NDOF + 1 : end);    
        AU   = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
        iL   = zeros(2 * KKT.NDOF, 1); 
        L    = zeros(2 * KKT.NDOF, 1);          
        return
    end        
end

if (fail) 
    fprintf('\nNewton method is converging, but did not reduce the error to the requested level\n')
    fprintf('\nIncrease number of Newton interations\n')   
    U    = V(1 : 2 * KKT.NDOF);
    Lamb = V(2 * KKT.NDOF + 1 : end);    
    AU   = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
    iL   = zeros(2 * KKT.NDOF, 1); 
    L    = zeros(2 * KKT.NDOF, 1);     
else
    if (flag == 0)
        AU = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
        iL = zeros(2 * KKT.NDOF, 1); 
        L  = zeros(2 * KKT.NDOF, 1);
    else
        [AU, iL] = computePartOfTimeDerivative1(msh, md, U, ts(m - 1) + dts(m - 1), flux, AE, A2, ME, IME, tc, quad1, bs);
        L = L + AU * U + iL; 
    end        
end

end