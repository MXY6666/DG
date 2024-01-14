% Compute the derivative of temperature with respect to the conserved 
% variables and w = \nabla rho
function dT = getTemperatureDeriv(rho, u, v, E, w1, w2, tc)

dT = zeros([size(rho), 6]);
dT(:, :, 1) = 3 / (8 * tc.Cv) * (3 - ((E - 0.5 / tc.We * (w1.^2 + w2.^2)) ./ rho - (u.^2 + v.^2)) ./ rho);
dT(:, :, 2) = -3 / (8 * tc.Cv) * u ./ rho;
dT(:, :, 3) = -3 / (8 * tc.Cv) * v ./ rho;
dT(:, :, 4) = 3 / (8 * tc.Cv) ./ rho;
dT(:, :, 5) = -3 / (8 * tc.Cv * tc.We) * w1 ./ rho;
dT(:, :, 6) = -3 / (8 * tc.Cv * tc.We) * w2 ./ rho;

end