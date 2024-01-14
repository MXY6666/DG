function [U, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, tc, quad, bs, bt, slo, KKT, newton, mtol)

fprintf('\n\nNewton solver for DIRK time integration\n')
% Initialize stage vector
KStage = zeros(3 * KKT.NDOF, bt.Knstages);
LStage = zeros(3 * KKT.NDOF, bt.Lnstages); 
if (bt.form == 1)
    KStage(:, 1) = U0;
end

stage.Ks = 1;
stage.Ls = 1;
for it = 0 : newton.restart
    % Set switch to adjust time step if convergence is very fast for all stages
    change_dt = false(bt.nstages, 1);

    % Solve each stage vector
    for s = 1 : bt.nstages
        % Store the stage index
        stage.s = s;    
        if (bt.form == 1)         
            stage.Ls = s;
        elseif (bt.form == 2)
            stage.Ks = s;                     
            if (s > 1)
                KStage(:, s) = KStage(:, s - 1);
            end
        end
        
        % Update stage vector
        if (KKT.enable)
            [KStage, Lamb, LStage, k, fail] = NewtonLeastSquaresSemiSmoothActiveSet(msh, md, U0, KStage, Lamb, LStage, nt, t, dt, stage, M, tc, bs, bt, KKT, newton, mtol);
        else
            [KStage, LStage, k, fail] = NewtonArmijoRule(msh, md, U0, KStage, LStage, nt, t, dt, stage, M, tc, bs, bt, KKT, newton, mtol);
        end
        
        % Judge if the Newton method converges well 
        if (fail)    
            if (it < newton.restart)
                fprintf('\nWARNING: Newton method is restarted due to lack of convergence\n')
                
                % Reduce time step
                dt = 0.5 * dt;
                
                % Initialize stage vector 
                if (bt.form == 1) 
                    KStage(:, 1) = U0;
                elseif (bt.form == 2)
                    KStage(:, 1) = 0;
                end
                
                break
            else
                error('\nERROR: Newton method was not able to converge after %d restarts\n', newton.restart)
            end
        else
            % Apply slope limiter
            KStage(:, stage.Ks) = applySlopeLimiter(msh, md, KStage(:, stage.Ks), tc, bs, slo);
            
            % Make sure the solution after slope limiter still satisfies
            % the constraints
            if (slo ~= 0 && KKT.enable)
                [KStage(:, 1), ~, fail] = constrainedProjection(msh, KStage(:, 1), nt, t + dt, 2, tc, quad, bs, KKT, newton, mtol);
                
                if (fail)
                    if (it < newton.restart)
                        fprintf('\nWARNING: Newton method is restarted due to lack of convergence in constrained projection\n')

                        % Reduce time step
                        dt = 0.75 * dt;

                        % Initialize stage vector 
                        KStage(:, 1) = U0;

                        break
                    else
                        error('\nERROR: Newton method was not able to converge after %d restarts\n', newton.restart)
                    end
                else
                    % Update the DG residual
                    if (s ~= bt.nstages) || ~(bt.type == 11 && any(bt.order == 1 : 4))
                        LStage(:, s) = computeDGResidual(msh, md, KStage(:, 1), tc, bs);
                    end
                end                
            end
                        
            % Increase time step if Newton method converges well
            if (k <= 10)
                change_dt(s) = true;
            end
        end  
    end
    
    if (~fail)
        % solution at new time level
        if (bt.form == 1)
            if (bt.type == 11 && any(bt.order == 1 : 4))
                U = KStage(:, 1);
            else
                U = U0 - dt * (M \ (LStage * bt.b)); 
            end
        elseif (bt.form == 2)
            U = U0 + dt * (KStage * bt.b);
        end
        
        if ~(bt.type == 11 && any(bt.order == 1 : 4))   
            % Apply slope limiter
            U = applySlopeLimiter(msh, md, U, tc, bs, slo);
            
            % Ensure that the DG solution at new time level satisfies the
            % constraints
            if (KKT.enable)
                [U, ~, fail] = constrainedProjection(msh, U, nt, t + dt, 2, tc, quad, bs, KKT, newton, mtol);
            end
            
            if (fail)
                if (it < newton.restart)
                    fprintf('\nWARNING: Newton method is restarted due to lack of convergence in constrained projection\n')

                    % Reduce time step
                    dt = 0.75 * dt;

                    % Initialize stage vector 
                    if (bt.form == 1) 
                        KStage(:, 1) = U0;
                    elseif (bt.form == 2)
                        KStage(:, 1) = 0;
                    end 
                    
                    continue
                else
                    error('\nERROR: Newton method was not able to converge after %d restarts\n', newton.restart)
                end
            end
        end    

        % new time level
        t = t + dt;

        % Compute global maximum and minimum 
        u = bs.phi{1} * reshape(U, [bs.nb, 3 * msh.nLElems]);
        p = getPressure(u(:, 1 : msh.nLElems), u(:, msh.nLElems + 1 : 2 * msh.nLElems), u(:, 2 * msh.nLElems + 1 : end), tc);
        rhomin = min(u(:, 1 : msh.nLElems), [], 'all');
        pmin   = min(p, [], 'all');
        fprintf('\n\nExtrema at time          %+15.6e\n', 0);
        fprintf('Minimum density          %+15.6e\n', rhomin);
        fprintf('Minimum pressure         %+15.6e\n', pmin);       
        break                       
    end
end

% Compute maximum time step at new time level
maxdt = setdt(msh, md, U, t, CFLMax, tp, tc, bs);

% Increase time step if Newton method converges well
if all(change_dt)
    dt = 1.2 * dt;
end
if (dt > maxdt)
    dt = maxdt;
end

end