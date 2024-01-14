% Plot solutions
function plotSols(msh, U, Q, Lamb, t, tc, quad, bs, KKT)

U = reshape(U, [bs.nb, msh.nLElems]);
Q = reshape(Q, [bs.nb, msh.nLElems]);
% mesh points
mp  = [msh.elemCenter(:, msh.LElems) - 0.5 * msh.elemLength(:, msh.LElems), tc.dm(2)];
mp  = sort(mp);

% Gauss quadrature points in each element
gp = msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad.points;

% exact solution at mesh points
u_mp = tc.u(mp, t);
q_mp = tc.u_x(mp, t);

% numerical solution at Gauss points
uh_gp = bs.phi{1} * U;
qh_gp = bs.phi{1} * Q;

% numerical solution at mesh points
if (quad.type == 101)
    F = griddedInterpolant(gp(:), uh_gp(:));
    uh_mp = F(mp);
    
    F = griddedInterpolant(gp(:), qh_gp(:));
    qh_mp = F(mp);
elseif (quad.type == 102)
    [~, I] = sort(gp(1, :));
    
    uh_gp = uh_gp(:, I);
    uh_mp = [uh_gp(1, 1), 0.5 * (uh_gp(end, 1 : end - 1) + uh_gp(1, 2 : end)), uh_gp(end, end)];
    
    qh_gp = qh_gp(:, I);
    qh_mp = [qh_gp(1, 1), 0.5 * (qh_gp(end, 1 : end - 1) + qh_gp(1, 2 : end)), qh_gp(end, end)];
end

% Plot
figure(1)
plot(mp, u_mp, 'b-', 'LineWidth', 2); hold on
plot(mp, uh_mp, 'b*');
legend('$u$', '$u_h$', 'Location','northeast', 'FontSize', 14, 'Interpreter', 'latex');
title(['t = ', num2str(t)], 'FontSize', 16)  
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$u/u_h$', 'FontSize', 20, 'Interpreter', 'latex');
hold off

figure(2)
plot(mp, q_mp, 'b-', 'LineWidth', 2); hold on
plot(mp, qh_mp, 'b*');
legend('$u_x$', '$q_h$', 'Location','southeast', 'FontSize', 14, 'Interpreter', 'latex');
title(['t = ', num2str(t)], 'FontSize', 16)  
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$u_x/q_h$', 'FontSize', 20, 'Interpreter', 'latex');
hold off

if (KKT.enable)
    inc = KKT.positivity + KKT.maximum;
    gp = repelem(gp, 1, inc);
    gp = gp(:);
    
    figure(3)
    plot(gp, Lamb, 'r*');
    title(['t = ', num2str(t)], 'FontSize', 16)
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    
    figure(4);
    colorL = [0 0 0];
    colorR = [0 0 0];
    set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   
    
    yyaxis 'left' 
    plot(mp, u_mp, 'b-', 'LineWidth', 2); hold on
    plot(mp, uh_mp, 'b*');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$u/u_h$', 'FontSize', 20, 'Interpreter', 'latex');
    axis([-0.5, 2, 0, 1.3])
    hold off
    
    pind  = Lamb > 1.e-10;
    pLamb = Lamb(pind);
    pgp   = gp(pind);
    
    yyaxis 'right'
    plot(pgp, pLamb, 'ro');
    ylabel('$\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    warning('off')
    legend('$u$', '$u_h$', '$\lambda$', 'Location','northeast', 'FontSize', 14, 'Interpreter', 'latex');
    warning('on')
end

end

