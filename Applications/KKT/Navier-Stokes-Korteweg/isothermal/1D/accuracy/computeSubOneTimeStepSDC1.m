% Compute one time step in a time subinterval 
function [U, AU, iL, L, it, fail] = computeSubOneTimeStepSDC1(msh, md, U0, AU0, iL0, L0, k, m, nt, ts, dts, flag, flux, AE, A2, A3, A32, ME, IME, tc, quad1, bs, SDC, KKT, newton)
    
NDOF = msh.nLElems * bs.nb;
t    = ts(1) + sum(dts);
% Initialize Newton solution vector
U = U0;

% Compute part of Jacobian matrix and residuals
if (k == 1)
    Jac1 = speye(2 * NDOF) - dts(m - 1) * AU0;
    res1 = -U0 - dts(m - 1) * iL0;    
else
    Jac1 = speye(2 * NDOF) - SDC.theta * dts(m - 1) * AU0;
    res1 = -U0 - SDC.theta * dts(m - 1) * (iL0 - L0(:, m)) - sum(dts) * sum(SDC.weights(m - 1, :) .* L0, 2); 
end

% Newton step
fail = true;
res = computeResidual(msh, md, U, [], AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
for it = 1 : newton.maxit
    if (it == 1 || L2_res > min(newton.tolres, 1.e-12))
        % Compute the Jacobian-associated matrix 
        Imp = computeImplicitMatrix(msh, md, U, [], [], [], Jac1, k, m, dts, [], flux, A3, A32, ME, IME, bs, SDC, KKT, newton);
        
        % Scale residual vectors
        scaRes = Imp.PR * res; 
        
        % Set switch to indicate that full Newton is used
        fullNewton = true;        
    else
        % Set switch to indicate that full Newton is used
        fullNewton = false;        
    end
    
    if (it == 1)        
        % Print convergence history 
        fprintf('\n\nTime-step      Time        k-m     Newton-iter   FN    L2-norm Residual   L2-norm DelU\n')  
        fprintf('%6d      %10.4e     %d-%d        %3d       %3d     %12.4e\n', nt, t, k, m - 1, 0, fullNewton, snorm(scaRes))
    end    
    
    % Solve linear system for Newton correction     
    scadU = linSolve(Imp.Mat, -scaRes, newton.linearSolver); 
    dU    = Imp.PC * scadU;
    
    % Store solution and residual vector of previous Newton iteration
    Uold      = U;
    scaResold = scaRes;
    
    % Compute trial solution
    U = Uold + dU;
    
    % Compute residual vector for new descent vector
    [res, iL] = computeResidual(msh, md, U, [], AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
    scaRes = Imp.PR * res;   
    
    % Compute L2-norm of residual
    L2_resold = snorm(scaResold);
    L2_res    = snorm(scaRes);
    
    % Perform line search based on Armijo rule
    if (L2_res >= (1 - newton.alpha) * L2_resold) 
        for l = 1 : newton.maxitLineSearch - 1
            % Compute trial solution
            U = Uold + newton.beta^l * dU;
            
            % Compute residual vector for new descent vector
            [res, iL] = computeResidual(msh, md, U, [], AU0, res1, k, m, dts, flux, A3, A32, ME, IME, bs, SDC, KKT);
            scaRes = Imp.PR * res;       
            
            % Compute L2-norm of residual
            L2_res = snorm(scaRes);
            
            if (L2_res < (1 - newton.alpha * newton.beta^l) * L2_resold)
                % Found proper search direction
                break
            end  
            
            if (l == newton.maxitLineSearch - 1)
                % Line search did not converge
                fprintf('\nNo feasible search direction found  -- restart Newton iteration with small time step\n')
                AU = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
                iL = zeros(2 * KKT.NDOF, 1); 
                L  = zeros(2 * KKT.NDOF, 1); 
                return
            end            
        end  
    end
    
    % Compute L2-norm of solution difference
    L2_dU = snorm(scadU);
    
    % Print convergence history 
    fprintf('%6d      %10.4e     %d-%d        %3d       %3d     %12.4e      %12.4e\n', nt, t, k, m - 1, it, fullNewton, L2_res, L2_dU)
    
    % Determine if stop the Newton method
    if (L2_res <= newton.tolres && L2_dU <= newton.toldV)
        fail = false;
        L = iL;
        break
    end
    
    % Restart Newton iteration if a local minimum is found
    if (L2_res > newton.tolres && L2_dU <= newton.toldVLocal)
        fprintf('\nLocal minimum found  -- restart Newton iteration with small time step\n')  
        AU = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
        iL = zeros(2 * KKT.NDOF, 1); 
        L  = zeros(2 * KKT.NDOF, 1);         
        return
    end        
end

if (fail) 
    fprintf('\nNewton method is converging, but did not reduce the error to the requested level\n')
    fprintf('\nIncrease number of Newton interations\n')
    AU = sparse([], [], [], 2 * KKT.NDOF, 2 * KKT.NDOF); 
    iL = zeros(2 * KKT.NDOF, 1); 
    L  = zeros(2 * KKT.NDOF, 1); 
else
    if (flag == 0)
        AU = sparse([], [], [], 2 * NDOF, 2 * NDOF); 
        iL = zeros(2 * NDOF, 1); 
        L  = zeros(2 * NDOF, 1);
    else
        [AU, iL] = computePartOfTimeDerivative1(msh, md, U, ts(m - 1) + dts(m - 1), flux, AE, A2, ME, IME, tc, quad1, bs);
        L = L + AU * U + iL; 
    end        
end

end