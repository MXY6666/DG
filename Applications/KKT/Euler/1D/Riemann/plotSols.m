% Plot solutions
function plotSols(msh, U, Lamb, t, pn, tc, quad, bs, KKT)

U = reshape(U, [bs.nb, 3 * msh.nLElems]);
% exact values to plot
if (tc.pro == 8)
    [xx, rho, u, p] = E1RPEX(1, 0.3, 1000, tc.gamma, t, tc.rhoL, tc.uL, tc.pL, tc.rhoR, tc.uR, tc.pR, 1);
elseif (tc.pro ~= 10 && tc.pro ~= 11)
    [rho, u, p] = computeExactSolution(pn.exa, t, tc);
end

% numerical values to plot
rhoh = 0.5 * reshape(quad.weights' * (bs.phi{1} * U(:, pn.numElems)), [length(pn.num), 1]);
mh   = 0.5 * reshape(quad.weights' * (bs.phi{1} * U(:, msh.nLElems + pn.numElems)), [length(pn.num), 1]);
Eh   = 0.5 * reshape(quad.weights' * (bs.phi{1} * U(:, 2 * msh.nLElems + pn.numElems)), [length(pn.num), 1]);
uh   = mh ./ rhoh;
ph   = getPressure(rhoh, mh, Eh, tc);

% Plot
if (KKT.positivity)
    gp = pn.exa(:);
    Lamb_rho = Lamb(1 : KKT.LDOF);
    Lamb_p   = Lamb(KKT.LDOF + 1 : end);
    
    pind      = Lamb_rho > 1.e-10;
    pgp_rho   = gp(pind);
    pLamb_rho = Lamb_rho(pind); 
    
    pind    = Lamb_p > 1.e-10;
    pgp_p   = gp(pind);
    pLamb_p = Lamb_p(pind);     
end
if (tc.pro == 9)
    xx = pn.exa;
end
if (tc.pro == 8 || tc.pro == 9)
    figure(1)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   

        yyaxis 'left' 
        semilogy(xx, rho, 'b-', 'LineWidth', 2); hold on
        semilogy(pn.num, rhoh, 'b*');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off

        yyaxis 'right'
        plot(pgp_rho, pLamb_rho, 'ro');
        ylabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$\rho$', '$\rho_h$', '$\lambda_{\rho}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on')        
    else   
        semilogy(xx, rho, 'b-', 'LineWidth', 2); hold on
        semilogy(pn.num, rhoh, 'b*');
        legend('$\rho$', '$\rho_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off
    end
    
    figure(2)
    plot(xx, u, 'b-', 'LineWidth', 2); hold on
    plot(pn.num, uh, 'b*');
    legend('$u$', '$u_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$u/u_h$', 'FontSize', 20, 'Interpreter', 'latex');
    hold off

    figure(3)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   

        yyaxis 'left' 
        semilogy(xx, p, 'b-', 'LineWidth', 2); hold on
        semilogy(pn.num, ph, 'b*');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off

        yyaxis 'right'
        plot(pgp_p, pLamb_p, 'ro');
        ylabel('$\lambda_{p}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$p$', '$p_h$', '$\lambda_{p}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on')  
    else     
        semilogy(xx, p, 'b-', 'LineWidth', 2); hold on
        semilogy(pn.num, ph, 'b*');
        legend('$p$', '$p_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off      
    end    
elseif (tc.pro == 10 || tc.pro == 11)
    figure(1)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);
        
        yyaxis 'left' 
        plot(pn.num, rhoh, 'b-*', 'LineWidth', 1.5);
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');

        yyaxis 'right'
        plot(pgp_rho, pLamb_rho, 'ro');
        ylabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$\rho_h$', '$\lambda_{\rho}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on') 
    else
        plot(pn.num, rhoh, 'b-*', 'LineWidth', 1.5);
        legend('$\rho_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
    end
    
    figure(2)
    plot(pn.num, uh, 'b-*', 'LineWidth', 1.5);
    legend('$u_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$u_h$', 'FontSize', 20, 'Interpreter', 'latex');

    figure(3)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);
        
        yyaxis 'left' 
        plot(pn.num, ph, 'b-*', 'LineWidth', 2);
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p_h$', 'FontSize', 20, 'Interpreter', 'latex');

        yyaxis 'right'
        plot(pgp_p, pLamb_p, 'ro');
        ylabel('$\lambda_{p}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$p$', '$p_h$', '$\lambda_{p}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on')   
    else
        plot(pn.num, ph, 'b-*', 'LineWidth', 2);
        legend('$p_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p_h$', 'FontSize', 20, 'Interpreter', 'latex');
    end
else
    figure(1)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   

        yyaxis 'left' 
        plot(pn.exa, rho, 'b-', 'LineWidth', 2); hold on
        plot(pn.num, rhoh, 'b*');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off

        yyaxis 'right'
        plot(pgp_rho, pLamb_rho, 'ro');
        ylabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$\rho$', '$\rho_h$', '$\lambda_{\rho}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on')        
    else   
        plot(pn.exa, rho, 'b-', 'LineWidth', 2); hold on
        plot(pn.num, rhoh, 'b*');
        legend('$\rho$', '$\rho_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$\rho/\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off
    end
    
    figure(2)
    plot(pn.exa, u, 'b-', 'LineWidth', 2); hold on
    plot(pn.num, uh, 'b*');
    legend('$u$', '$u_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
    title(['t = ', num2str(t)], 'FontSize', 16)  
    xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$u/u_h$', 'FontSize', 20, 'Interpreter', 'latex');
    hold off

    figure(3)
    if (KKT.positivity)
        colorL = [0 0 0];
        colorR = [0 0 0];
        set(gcf, 'defaultAxesColorOrder', [colorL; colorR]);   

        yyaxis 'left' 
        plot(pn.exa, p, 'b-', 'LineWidth', 2); hold on
        plot(pn.num, ph, 'b*');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off

        yyaxis 'right'
        plot(pgp_p, pLamb_p, 'ro');
        ylabel('$\lambda_{p}$', 'FontSize', 20, 'Interpreter', 'latex');
        warning('off')
        legend('$p$', '$p_h$', '$\lambda_{p}$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        warning('on')  
    else     
        plot(pn.exa, p, 'b-', 'LineWidth', 2); hold on
        plot(pn.num, ph, 'b*');
        legend('$p$', '$p_h$', 'Location','northwest', 'FontSize', 14, 'Interpreter', 'latex');
        title(['t = ', num2str(t)], 'FontSize', 16)  
        xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
        ylabel('$p/p_h$', 'FontSize', 20, 'Interpreter', 'latex');
        hold off      
    end
end

end