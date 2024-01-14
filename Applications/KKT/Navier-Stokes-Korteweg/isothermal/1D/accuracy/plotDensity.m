% Plot solutions
function plotDensity(msh, U, Lamb, t, tc, quad1, bs, KKT)

U = reshape(U, [bs.nb, 2 * msh.nLElems]);
% mesh points
mp = [msh.elemCenter(:, msh.LElems) - 0.5 * msh.elemLength(:, msh.LElems), tc.dm(2)];
mp = sort(mp);

% Gauss quadrature points in each element
gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points;

% exact density at mesh points
rho_mp = tc.rho(mp, t);

% numerical density at Gauss points
rhoh_gp = bs.phi{1} * U(:, 1 : msh.nLElems);

% numerical density at mesh points
if (quad1.type == 101)
    F = griddedInterpolant(gp(:), rhoh_gp(:));
    rhoh_mp = F(mp);
elseif (quad1.type == 102)
    [~, I] = sort(gp(1, :));
    rhoh_gp = rhoh_gp(:, I);
    rhoh_mp = [0.5 * (rhoh_gp(1, 1) + rhoh_gp(end, end)), 0.5 * (rhoh_gp(end, 1 : end - 1) + rhoh_gp(1, 2 : end)), 0.5 * (rhoh_gp(1, 1) + rhoh_gp(end, end))];
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

if (KKT.positivity)
    gp = gp(:);
    figure(2)
    plot(gp, Lamb, 'r*');
    title(['t = ', num2str(t)], 'FontSize', 16)
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    
    figure(3);
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
    
    pind  = Lamb > 1.e-10;
    pLamb = Lamb(pind);
    pgp   = gp(pind);
    
    yyaxis 'right'
    plot(pgp, pLamb, 'ro');
    ylabel('$\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    warning('off')
    legend('$\rho$', '$\rho_h$', '$\lambda$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    warning('on')
end

end

