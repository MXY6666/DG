function [U, Lamb, t, dt] = computeOneTimeStepSDC(msh, md, U0, Lamb, nt, t, dt, maxdt, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton)

if (dt > maxdt)
    dt = maxdt;
end
fprintf('\n\nNewton solver for SDC time integration\n')
% initialization
NDOF = msh.nLElems * bs.nb;
AU = cell(1, SDC.order);
iL = zeros(2 * NDOF, SDC.order);

% values at first point
[AU{1}, iL(:, 1)] = computePartOfTimeDerivative1(msh, md, U0, t, flux, AE, A2, ME, IME, tc, quad1, bs);
if (SDC.order > 1)
    [~, iL0]  = computePartOfTimeDerivative2(msh, md, U0, 2, flux, ME, IME, bs);
    L0 = zeros(2 * NDOF, SDC.order);
    L  = zeros(2 * NDOF, SDC.order);
    L0(:, 1) = AU{1} * U0 + iL(:, 1) + iL0; 
    L( :, 1) = L0(:, 1);
end

for it = 0 : newton.restart
    % Set switch to adjust time step if convergence is very fast 
    if (SDC.order == 1)
        change_dt = false;
    else
        change_dt = false((SDC.order - 1) * SDC.order, 1);
        cnt = 0;
    end

    % times and time subintervals
    ts  = t + 0.5 * dt * (1 + SDC.points);
    dts = diff(ts);    
    
    if (SDC.order == 1)
        if (KKT.enable)
            [U, Lamb, ~, ~, ~, nit, fail] = computeSubOneTimeStepSDC2(msh, md, U0, Lamb, AU{1}, iL(:, 1), [], 1, 2, nt, ts, dts, 0, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);
        else
            [U, ~, ~, ~, nit, fail] = computeSubOneTimeStepSDC1(msh, md, U0, AU{1}, iL(:, 1), [], 1, 2, nt, ts, dts, 0, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);
        end
        
        % Judge if the Newton method converges well 
        if (fail)
            fprintf('\nWARNING: Newton method is restarted due to lack of convergence\n')

            % Reduce time step
            dt = 0.5 * dt;
            break
        else
            % Increase time step if Newton method converges well
            if (nit <= 10 && 1.2 * dt <= maxdt)
                change_dt = true;
            end
        end        
    else
        % first order approximation
        for m = 2 : SDC.order
            if (m == 2)
                U = U0;
            end
            if (KKT.enable)
                [U, Lamb, AU{m}, iL(:, m), L0(:, m), nit, fail] = computeSubOneTimeStepSDC2(msh, md, U, Lamb, AU{m - 1}, iL(:, m - 1), [], 1, m, nt, ts, dts, 1, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);
            else
                [U, AU{m}, iL(:, m), L0(:, m), nit, fail] = computeSubOneTimeStepSDC1(msh, md, U, AU{m - 1}, iL(:, m - 1), [], 1, m, nt, ts, dts, 1, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);  
            end
            
            % Judge if the Newton method converges well 
            if (fail)
                fprintf('\nWARNING: Newton method is restarted due to lack of convergence\n')

                % Reduce time step
                dt = 0.5 * dt;
                break
            else
                % Increase time step if Newton method converges well
                if (nit <= 10 && 1.2 * dt <= maxdt)
                    cnt = cnt + 1;
                    change_dt(cnt) = true;
                end
            end            
        end 

        % successive corrections
        if (~fail)
            for k = 2 : SDC.order
                if (k < SDC.order)
                    flag = 1;
                else
                    flag = 0;
                end

                for m = 2 : SDC.order
                    if (m == 2)
                        U = U0;
                    end
                    if (KKT.enable)
                        [U, Lamb, AU{m}, iL(:, m), L(:, m), nit, fail] = computeSubOneTimeStepSDC2(msh, md, U, Lamb, AU{m}, iL(:, m), L0, k, m, nt, ts, dts, flag, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton); 
                    else
                        [U, AU{m}, iL(:, m), L(:, m), nit, fail] = computeSubOneTimeStepSDC1(msh, md, U, AU{m}, iL(:, m), L0, k, m, nt, ts, dts, flag, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);
                    end
                    
                    % Judge if the Newton method converges well 
                    if (fail)
                        fprintf('\nWARNING: Newton method is restarted due to lack of convergence\n')

                        % Reduce time step
                        dt = 0.5 * dt;
                        break
                    else
                        % Increase time step if Newton method converges well
                        if (nit <= 10 && 1.2 * dt <= maxdt)
                            cnt = cnt + 1;
                            change_dt(cnt) = true;
                        end
                    end                     
                end
                
                if (fail)
                    break
                else
                    if (k < SDC.order)
                        L0 = L;
                    end
                end
            end
        end
    end
    
    if (~fail)
        % new time level
        t = t + dt;   
        
        % Compute global maximum density
        rho = bs.phi{1} * reshape(U(1 : KKT.NDOF), [bs.nb, msh.nLElems]);
        rhomin = min(rho, [], 'all');
        fprintf('\n\nExtrema at time          %+15.6e\n', t)
        fprintf('Minimum density          %+15.6e\n', rhomin)        
        break
    end
end

% Increase time step if Newton method converges well
if all(change_dt)
    dt = 1.2 * dt;
end

end

