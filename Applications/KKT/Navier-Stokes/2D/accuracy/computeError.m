function [L2err, Lierr] = computeError(msh, RHO, M, E, t, tc, quad2, bs)

L2err = zeros(3, 1);
Lierr = zeros(3, 1);
% Eval exact solutions at Gauss quadrature points
ct  = msh.elemCenter(:, msh.LElems);
h   = msh.elemLength(:, msh.LElems);
gpx = ct(1, :) + h(1, :) / 2 .* quad2.points(:, 1);
gpy = ct(2, :) + h(2, :) / 2 .* quad2.points(:, 2);
rho = tc.rho(gpx, gpy, t);
u   = [tc.u(gpx, gpy, t), tc.v(gpx, gpy, t)];
p   = tc.p(gpx, gpy, t);

% Eval numerical solutions at Gauss quadrature points
rhoh = bs.phi{1} * RHO;
mh   = bs.phi{1} * M;
Eh   = bs.phi{1} * E;
uh   = mh ./ repmat(rhoh, [1, 2]);
ph   = getPressure(rhoh, mh(:, 1 : msh.nLElems), mh(:, msh.nLElems + 1 : end), Eh, tc);

% compute L2 error
J = msh.elemJac(:, msh.LElems);
L2err(1) = sqrt(sum(J .* (quad2.weights' * (rho - rhoh).^2)));
L2err(2) = sqrt(sum(repmat(J, [1, 2]) .* (quad2.weights' * (u - uh).^2)));
L2err(3) = sqrt(sum(J .* (quad2.weights' * (p - ph).^2)));

% compute L_infinity error
Lierr(1) = max(max(abs(rho - rhoh)));  
Lierr(2) = max(max(abs(u - uh))); 
Lierr(3) = max(max(abs(p - ph))); 

end