% Compute the derivative of pressure with respect to the conserved 
% variables and w = \nabla rho
function dp = getPressureDeriv(rho, u, v, E, w1, w2, tc)

dp = zeros([size(rho), 6]);
dp(:, :, 1) = 3 / tc.Cv ./ (3 - rho).^2 .* (3 * rho .* (6 - rho) + (1.5 - rho) .* (u.^2 + v.^2) + E - 0.5 / tc.We * (w1.^2 + w2.^2)) - 6 * rho; 
dp(:, :, 2) = 3 / tc.Cv ./ (rho - 3) .* u;
dp(:, :, 3) = 3 / tc.Cv ./ (rho - 3) .* v;
dp(:, :, 4) = 3 / tc.Cv ./ (3 - rho);
dp(:, :, 5) = 3 / tc.Cv / tc.We ./ (rho - 3) .* w1;
dp(:, :, 6) = 3 / tc.Cv / tc.We ./ (rho - 3) .* w2;

end