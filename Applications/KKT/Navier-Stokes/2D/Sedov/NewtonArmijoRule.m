% Solve one stage with Newton technique for DIRK time integration method
function [KStage, LStage, k, fail] = NewtonArmijoRule(msh, md, U0, KStage, LStage, nt, t, dt, stage, isViscous, M, pr12, pc12, pv12, IME, tc, bs, bt, KKT, newton, mtol)

% Initialize Newton solution vector
U = KStage(:, stage.Ks); 

% Newton step
[res, ~, waveVel] = computeResidual(msh, md, U0, KStage, [], LStage, dt, stage, isViscous, M, IME, tc, bs, bt, KKT);
fail = true;
for k = 1 : newton.maxit
    if (k == 1 || L2_res > min(newton.tolres, 1.e-12))
        % Compute the Jacobian matrix 
        Imp = computeImplicitMatrix(msh, md, U0, KStage, [], [], [], dt, stage, [], isViscous, pr12, pc12, pv12, IME, tc, bs, bt, KKT, newton, mtol);

        % Scale residual vectors
        scaRes = Imp.PR * res;

        % Set switch to indicate that full Newton is used
        fullNewton = true;   
    else
        % Set switch to indicate that full Newton is used
        fullNewton = false; 
    end

    if (k == 1)
        % Compute CFL number
        CFL = waveVel * dt / md.meshSize;
        
        % Print convergence history 
        fprintf('\n\nTime-step      Time       RK-stage     Newton-iter      CFL       FN    L2-norm Residual   L2-norm DelV\n')  
        fprintf('%6d      %10.4e      %3d           %3d       %10.4e   %3d      %12.4e\n', nt, t + dt, stage.s, 0, CFL, fullNewton, snorm(scaRes))       
    end        

    % Solve linear system for Newton correction
    if (newton.linearSolver.type == 1)
        scadU = Imp.Mat \ (-scaRes);
    else
        scadU = linSolve(Imp.Mat, -scaRes, newton.linearSolver);
    end
    dU = Imp.PC * scadU;
    
    % Store solution and residual vector of previous Newton iteration
    Uold      = U;
    scaResold = scaRes;    

    % Compute trial solution
    U = Uold + dU;
    KStage(:, stage.Ks) = U;

    % Compute residual vector for new descent vector
    [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, [], LStage, dt, stage, isViscous, M, IME, tc, bs, bt, KKT);
    scaRes = Imp.PR * res;
    
    % Compute L2-norm of residual
    L2_resold = snorm(scaResold);
    L2_res    = snorm(scaRes);
    
    % Perform line search based on Armijo rule
    if (L2_res >= (1 - newton.alpha) * L2_resold) 
        for l = 1 : newton.maxitLineSearch - 1
            % Compute trial solution
            U = Uold + newton.beta^l * dU;
            KStage(:, stage.Ks) = U;
            
            % Compute residual vector for new descent vector
            [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, [], LStage, dt, stage, isViscous, M, IME, tc, bs, bt, KKT);
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
                return
            end            
        end  
    end
    
    % Compute CFL number
    CFL = waveVel * dt / md.meshSize;
 
    % Compute L2-norm of solution difference
    L2_dU = snorm(scadU);

    % Print convergence history 
    fprintf('%6d      %10.4e      %3d           %3d       %10.4e   %3d      %12.4e    %12.4e\n', nt, t + dt, stage.s, k, CFL, fullNewton, L2_res, L2_dU) 

    % Determine if stop the Newton method
    if (L2_res <= newton.tolres && L2_dU <= newton.toldV)
        fail = false;
        return
    end

    % Restart Newton iteration if a local minimum is found
    if (L2_res > newton.tolres && L2_dU <= newton.toldVLocal)
        fprintf('\nLocal minimum found  -- restart Newton iteration with small time step\n')  
        return
    end        
end    
fprintf('\nNewton method is converging, but did not reduce the error to the requested level\n')
fprintf('\nIncrease number of Newton interations\n')
  
end
