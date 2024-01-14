% Plot solutions
function plotSols(msh, U, Lamb, t, tc, quad2, bs, KKT)

U = reshape(U, [bs.nb, 4 * msh.nLElems]);
% mesh points
mpx = repmat(linspace(msh.dm(1), msh.dm(2), msh.N(1) + 1)', [1, msh.N(2) + 1]);
mpy = repmat(linspace(msh.dm(3), msh.dm(4), msh.N(2) + 1), [msh.N(1) + 1, 1]);

% Gauss quadrature points in each element
gpx = msh.elemCenter(1, msh.LElems) + 0.5 * msh.elemLength(1, msh.LElems) .* quad2.points(:, 1);
gpy = msh.elemCenter(2, msh.LElems) + 0.5 * msh.elemLength(2, msh.LElems) .* quad2.points(:, 2);

% numerical solution at Gauss points
uh_gp   = bs.phi{1} * U;
rhoh_gp = uh_gp(:, 1 : msh.nLElems);
ph_gp   = getPressure(uh_gp(:, 1 : msh.nLElems), uh_gp(:, msh.nLElems + 1 : 2 * msh.nLElems), uh_gp(:, 2 * msh.nLElems + 1 : 3 * msh.nLElems), uh_gp(:, 3 * msh.nLElems + 1 : end), tc);

% numerical solution at mesh points
rhoh_mp = griddata(gpx, gpy, rhoh_gp, mpx, mpy);
ph_mp   = griddata(gpx, gpy, ph_gp, mpx, mpy);

% Plot numerical solution
figure(1)
surf(mpx, mpy, rhoh_mp); 
title(['t = ', num2str(t)], 'FontSize', 16);  
xlabel('$x$',      'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$',      'FontSize', 20, 'Interpreter', 'latex');
zlabel('$\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

figure(2)
surf(mpx, mpy, ph_mp); 
title(['t = ', num2str(t)], 'FontSize', 16);  
xlabel('$x$',   'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$',   'FontSize', 20, 'Interpreter', 'latex');
zlabel('$p_h$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

if (KKT.positivity)
    gpx = gpx(:);
    gpy = gpy(:);
    Lamb_rho = Lamb(1 : KKT.LDOF);
    Lamb_p   = Lamb(KKT.LDOF + 1 : end);
    
    figure(3);
    pind  = Lamb_rho > 1.e-10;
    pLamb = Lamb_rho(pind);
    pgpx  = gpx(pind);
    pgpy  = gpy(pind);
    
    stem3(pgpx, pgpy, pLamb, ':*r');
    title(['t = ', num2str(t)], 'FontSize', 16);  
    xlabel('$x$',              'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$y$',              'FontSize', 20, 'Interpreter', 'latex');
    zlabel('$\lambda_{\rho}$', 'FontSize', 20, 'Interpreter', 'latex');
    axis(msh.dm)
    az = 55;
    el = 35;
    view(az, el);    
    
    figure(4);
    pind  = Lamb_p > 1.e-10;
    pLamb = Lamb_p(pind);
    pgpx  = gpx(pind);
    pgpy  = gpy(pind);
    
    stem3(pgpx, pgpy, pLamb, ':*r');
    title(['t = ', num2str(t)], 'FontSize', 16);  
    xlabel('$x$',         'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$y$',         'FontSize', 20, 'Interpreter', 'latex');
    zlabel('$\lambda_p$', 'FontSize', 20, 'Interpreter', 'latex');
    axis(msh.dm)
    az = 55;
    el = 35;
    view(az, el);
end

end

