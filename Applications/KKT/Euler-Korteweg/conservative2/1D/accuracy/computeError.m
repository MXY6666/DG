function [L2err, Lierr] = computeError(msh, U, t, tc, quad1, bs)

L2err = zeros(2, 1);
Lierr = zeros(2, 1);
% Eval exact solutions at Gauss quadrature points
ct  = msh.elemCenter(:, msh.LElems);
h   = msh.elemLength(:, msh.LElems);
gpx = ct + 0.5 * h .* quad1.points;
rho = tc.rho(gpx, t);
u   = tc.u(gpx, t);

% Eval numerical solutions at Gauss quadrature points
rhoh = bs.phi{1} * U(:, 1 : msh.nLElems);
mh   = bs.phi{1} * U(:, msh.nLElems + 1 : end);
uh   = mh ./ rhoh;

% Compute L2 error
J = msh.elemJac(:, msh.LElems);
L2err(1) = sqrt(sum(J .* (quad1.weights' * (rho - rhoh).^2)));
L2err(2) = sqrt(sum(J .* (quad1.weights' * (u - uh).^2)));

% Compute L-infinity error
Lierr(1) = max(max(abs(rho - rhoh)));  
Lierr(2) = max(max(abs(u - uh))); 

end