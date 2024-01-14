% Plot solutions
function plotSols(msh, U, Lamb, t, tc, quad2, bs, KKT)

U = reshape(U, [bs.nb, msh.nLElems]);
% mesh points
mpx = repmat(linspace(msh.dm(1), msh.dm(2), msh.N(1) + 1)', [1, msh.N(2) + 1]);
mpy = repmat(linspace(msh.dm(3), msh.dm(4), msh.N(2) + 1), [msh.N(1) + 1, 1]);

% Gauss quadrature points in each element
gpx = msh.elemCenter(1, msh.LElems) + 0.5 * msh.elemLength(1, msh.LElems) .* quad2.points(:, 1);
gpy = msh.elemCenter(2, msh.LElems) + 0.5 * msh.elemLength(2, msh.LElems) .* quad2.points(:, 2);

% exact solution at mesh points
u_mp = tc.u(mpx, mpy, t);

% numerical solution at Gauss points
uh_gp = bs.phi{1} * U;

% numerical solution at mesh points
uh_mp = griddata(gpx, gpy, uh_gp, mpx, mpy);

% Plot analytic solution
figure(1)
surf(mpx, mpy, u_mp); 
title(['t = ', num2str(t)], 'FontSize', 16); 
xlabel('$x$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 20, 'Interpreter', 'latex');
zlabel('$u$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

% Plot numerical solution
figure(2)
surf(mpx, mpy, uh_mp); 
title(['t = ', num2str(t)], 'FontSize', 16);  
xlabel('$x$',   'FontSize', 20, 'Interpreter', 'latex');
ylabel('$y$',   'FontSize', 20, 'Interpreter', 'latex');
zlabel('$u_h$', 'FontSize', 20, 'Interpreter', 'latex');
az = 55;
el = 35;
view(az, el);

if (KKT.enable)
    inc = KKT.positivity + KKT.maximum;
    gpx = repelem(gpx, 1, inc);
    gpy = repelem(gpy, 1, inc);
    gpx = gpx(:);
    gpy = gpy(:);
    
    figure(3);
    pind  = Lamb > 1.e-10;
    pLamb = Lamb(pind);
    pgpx  = gpx(pind);
    pgpy  = gpy(pind);
    
    stem3(pgpx, pgpy, pLamb, ':*r');
    title(['t = ', num2str(t)], 'FontSize', 16);  
    xlabel('$x$',       'FontSize', 20, 'Interpreter', 'latex');
    ylabel('$y$',       'FontSize', 20, 'Interpreter', 'latex');
    zlabel('$\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    axis(msh.dm)
    az = 55;
    el = 35;
    view(az, el);    
end

end

