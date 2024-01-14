% Plot solutions
function plotSols(msh, U, W1, W2, Lamb, t, tc, quad2, bs, KKT)

U  = reshape(U, [bs.nb, 4 * msh.nLElems]);
W1 = reshape(W1, [bs.nb, msh.nLElems]);
W2 = reshape(W2, [bs.nb, msh.nLElems]);
% mesh points
mpx = repmat(linspace(msh.dm(1), msh.dm(2), msh.N(1) + 1)', [1, msh.N(2) + 1]);
mpy = repmat(linspace(msh.dm(3), msh.dm(4), msh.N(2) + 1), [msh.N(1) + 1, 1]);

% Gauss quadrature points in each element
gpx = msh.elemCenter(1, msh.LElems) + 0.5 * msh.elemLength(1, msh.LElems) .* quad2.points(:, 1);
gpy = msh.elemCenter(2, msh.LElems) + 0.5 * msh.elemLength(2, msh.LElems) .* quad2.points(:, 2);

% exact solution at mesh points
rho_mp = tc.rho(mpx, mpy, t);
T_mp   = tc.T(mpx, mpy, t);
p_mp   = getPressure(rho_mp, T_mp);

% numerical solution at Gauss points
uh_gp   = bs.phi{1} * U;
w1h_gp  = bs.phi{1} * W1;
w2h_gp  = bs.phi{1} * W2;
rhoh_gp = uh_gp(:, 1 : msh.nLElems);
vh_gp   = uh_gp(:, msh.nLElems + 1 : 3 * msh.nLElems) ./ repmat(rhoh_gp, [1, 2]);
Th_gp   = getTemperature(uh_gp(:, 1 : msh.nLElems), vh_gp(:, 1 : end / 2), vh_gp(:, end / 2 + 1 : end), uh_gp(:, 3 * msh.nLElems + 1 : end), w1h_gp, w2h_gp, tc);
ph_gp   = getPressure(rhoh_gp, Th_gp);

% numerical solution at mesh points
rhoh_mp = griddata(gpx, gpy, rhoh_gp, mpx, mpy);
ph_mp   = griddata(gpx, gpy, ph_gp, mpx, mpy);

% Plot analytic solution
figure(1)
surf(mpx, mpy, rho_mp); 
title(['t = ', num2str(t)], 'FontSize', 16); 
xlabel('$x$',    'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$',    'FontSize', 20, 'Interpreter', 'latex');
zlabel('$\rho$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

figure(2)
surf(mpx, mpy, p_mp); 
title(['t = ', num2str(t)], 'FontSize', 16); 
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 20, 'Interpreter', 'latex');
zlabel('$p$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

% Plot numerical solution
figure(3)
surf(mpx, mpy, rhoh_mp); 
title(['t = ', num2str(t)], 'FontSize', 16);  
xlabel('$x$',      'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$',      'FontSize', 20, 'Interpreter', 'latex');
zlabel('$\rho_h$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

figure(4)
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
    
    figure(5);
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
    
    figure(6);
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

