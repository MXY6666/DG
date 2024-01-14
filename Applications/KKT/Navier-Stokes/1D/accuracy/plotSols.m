% Plot solutions
function plotSols(msh, U, Lamb, t, tc, quad1, bs, KKT)

U = reshape(U, [bs.nb, 3 * msh.nLElems]);
% mesh points
mp = [msh.elemCenter(:, msh.LElems) - 0.5 * msh.elemLength(:, msh.LElems), tc.dm(2)];
mp = sort(mp);

% Gauss quadrature points in each element
gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points;

% exact solutions at mesh points
rho_mp = tc.rho(mp, t);
p_mp   = tc.p(mp, t);

% numerical solutions at Gauss points
uh_gp   = bs.phi{1} * U;
rhoh_gp = uh_gp(:, 1 : msh.nLElems);
ph_gp   = getPressure(uh_gp(:, 1 : msh.nLElems), uh_gp(:, msh.nLElems + 1 : 2 * msh.nLElems), uh_gp(:, 2 * msh.nLElems + 1 : end), tc);

% numerical solutions at mesh points
if (quad1.type == 101)
    F_rho = griddedInterpolant(gp(:), rhoh_gp(:));
    F_p   = griddedInterpolant(gp(:), ph_gp(:));
    rhoh_mp = F_rho(mp);
    ph_mp   = F_p(mp);
elseif (quad1.type == 102)
    [~, I] = sort(gp(1, :));
    
    rhoh_gp = rhoh_gp(:, I);
    ph_gp   = ph_gp(:, I);
    rhoh_mp = [0.5 * (rhoh_gp(1, 1) + rhoh_gp(end, end)), 0.5 * (rhoh_gp(end, 1 : end - 1) + rhoh_gp(1, 2 : end)), 0.5 * (rhoh_gp(1, 1) + rhoh_gp(end, end))];
    ph_mp   = [0.5 * (ph_gp(1, 1) + ph_gp(end, end)), 0.5 * (ph_gp(end, 1 : end - 1) + ph_gp(1, 2 : end)), 0.5 * (ph_gp(1, 1) + ph_gp(end, end))];    
end

% Plot
figure(1)
plot(mp, rho_mp, 'b-', 'LineWidth', 2); hold on
plot(mp, rhoh_mp, 'b*');
legend('$\rho$', '$\rho_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
title(['t = ', num2str(t)], 'FontSize', 16)  
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
hold off

figure(2)
plot(mp, p_mp, 'b-', 'LineWidth', 2); hold on
plot(mp, ph_mp, 'b*');
legend('$p$', '$p_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
title(['t = ', num2str(t)], 'FontSize', 16)  
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
hold off

if (KKT.positivity)
    gp = gp(:);
    Lamb_rho = Lamb(1 : KKT.LDOF);
    Lamb_p   = Lamb(KKT.LDOF + 1 : end);
    
    figure(3)
    plot(gp, Lamb_rho, 'r*');
    title(['t = ', num2str(t)], 'FontSize', 16)
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
    
    figure(4)
    plot(gp, Lamb_p, 'r*');
    title(['t = ', num2str(t)], 'FontSize', 16)
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$\lambda_{p}$', 'FontSize', 20, 'Interpreter', 'latex');    
    
    figure(5);
    colorL = [0 0 0];
    colorR = [0 0 0];
    set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   
    
    yyaxis 'left' 
    plot(mp, rho_mp, 'b-', 'LineWidth', 2); hold on
    plot(mp, rhoh_mp, 'b*');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
    hold off
    
    pind  = Lamb_rho > 1.e-10;
    pgp   = gp(pind);
    pLamb = Lamb_rho(pind); 
    
    yyaxis 'right'
    plot(pgp, pLamb, 'ro');
    ylabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
    warning('off')
    legend('$\rho$', '$\rho_h$', '$\lambda_{\rho}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    warning('on')
    
    figure(6);
    colorL = [0 0 0];
    colorR = [0 0 0];
    set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   
    
    yyaxis 'left' 
    plot(mp, p_mp, 'b-', 'LineWidth', 2); hold on
    plot(mp, ph_mp, 'b*');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
    hold off
    
    pind  = Lamb_p > 1.e-10;
    pgp   = gp(pind);
    pLamb = Lamb_p(pind); 
    
    yyaxis 'right'
    plot(pgp, pLamb, 'ro');
    ylabel('$\lambda_{p}$', 'FontSize', 20, 'Interpreter', 'latex');
    warning('off')
    legend('$p$', '$p_h$', '$\lambda_{p}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    warning('on')    
end

end
