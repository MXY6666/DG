function [U, W1, W2, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, Ax, Ay, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, quad2, bs, bt, inviscidFlux, KKT, newton, mtol)

fprintf('\n\nNewton solver for DIRK time integration\n')
% Initialize stage vector
KStage = zeros(4 * KKT.NDOF, bt.Knstages);
LStage = zeros(4 * KKT.NDOF, bt.Lnstages); 
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
            [KStage, Lamb, LStage, k, fail] = NewtonLeastSquaresSemiSmoothActiveSet(msh, md, U0, KStage, Lamb, LStage, nt, t, dt, stage, M, Ax, Ay, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, quad2, bs, bt, inviscidFlux, KKT, newton, mtol);
        else
            [KStage, LStage, k, fail] = NewtonArmijoRule(msh, md, U0, KStage, LStage, nt, t, dt, stage, M, Ax, Ay, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, quad2, bs, bt, inviscidFlux, KKT, newton, mtol);
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
        
        % Ensure that the DG solution at new time level satisfies the
        % constraints
        if (KKT.enable)
            if ~(bt.type == 11 && any(bt.order == 1 : 4))
                [U, ~, fail] = constrainedProjection(msh, U, nt, t + dt, 2, Ax, Ay, tc, quad2, bs, KKT, newton, mtol);
            end
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
            else
                error('\nERROR: Newton method was not able to converge after %d restarts\n', newton.restart)
            end            
        else
            % new time level
            t = t + dt;
            
            % Compute global maximum and minimum 
            W1 = Ax * U(1 : end / 4); 
            W2 = Ay * U(1 : end / 4);
            u  = bs.phi{1} * reshape(U, [bs.nb, 4 * msh.nLElems]);
            w1 = bs.phi{1} * reshape(W1, [bs.nb, msh.nLElems]);
            w2 = bs.phi{1} * reshape(W2, [bs.nb, msh.nLElems]);
            v  = u(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(u(:, 1 : msh.nLElems), [1, 2]);
            T = getTemperature(u(:, 1 : msh.nLElems), v(:, 1 : end / 2), v(:, end / 2 + 1 : end), u(:, 3 * msh.nLElems + 1 : end), w1, w2, tc);
            p = getPressure(u(:, 1 : msh.nLElems), T);
            rhomin = min(u(:, 1 : msh.nLElems), [], 'all');
            pmin   = min(p, [], 'all');
            fprintf('\n\nExtrema at time          %+15.6e\n', 0);
            fprintf('Minimum density          %+15.6e\n', rhomin);
            fprintf('Minimum pressure         %+15.6e\n', pmin);       
            break            
        end              
    end
end

% Compute maximum time step at new time level
maxdt = setdt(msh, md, U, W1, W2, t, CFLMax, tp, tc, bs);

% Increase time step if Newton method converges well
if all(change_dt)
    dt = 1.2 * dt;
end
if (dt > maxdt)
    dt = maxdt;
end

end
