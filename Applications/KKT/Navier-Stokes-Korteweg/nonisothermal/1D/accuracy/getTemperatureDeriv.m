% Compute the derivative of temperature with respect to the conserved 
% variables and w = \nabla rho
function dT = getTemperatureDeriv(rho, u, E, w, tc)

dT = zeros([size(rho), 4]);
dT(:, :, 1) = 3 / (8 * tc.Cv) * (3 - ((E - 0.5 / tc.We * w.^2) ./ rho - u.^2) ./ rho);
dT(:, :, 2) = -3 / (8 * tc.Cv) * u ./ rho;
dT(:, :, 3) = 3 / (8 * tc.Cv) ./ rho;
dT(:, :, 4) = -3 / (8 * tc.Cv * tc.We) * w ./ rho;

end