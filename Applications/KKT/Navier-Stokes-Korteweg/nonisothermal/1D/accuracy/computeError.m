function [L2err, Lierr] = computeError(msh, RHO, M, E, W, t, tc, quad1, bs)

L2err = zeros(3, 1);
Lierr = zeros(3, 1);
% Eval exact solutions at Gauss quadrature points
ct  = msh.elemCenter(:, msh.LElems);
h   = msh.elemLength(:, msh.LElems);
gpx = ct + 0.5 * h .* quad1.points;
rho = tc.rho(gpx, t);
u   = tc.u(gpx, t);
p   = tc.p(gpx, t);

% Eval numerical solutions at Gauss quadrature points
rhoh = bs.phi{1} * RHO;
mh   = bs.phi{1} * M;
Eh   = bs.phi{1} * E;
wh   = bs.phi{1} * W;
uh   = mh ./ rhoh;
Th   = getTemperature(rhoh, mh ./ rhoh, Eh, wh, tc);
ph   = getPressure(rhoh, Th);

% Compute L2 error
J = msh.elemJac(:, msh.LElems);
L2err(1) = sqrt(sum(J .* (quad1.weights' * (rho - rhoh).^2)));
L2err(2) = sqrt(sum(J .* (quad1.weights' * (u - uh).^2)));
L2err(3) = sqrt(sum(J .* (quad1.weights' * (p - ph).^2)));

% Compute L-infinity error
Lierr(1) = max(max(abs(rho - rhoh)));  
Lierr(2) = max(max(abs(u - uh))); 
Lierr(3) = max(max(abs(p - ph))); 

 
end