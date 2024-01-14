% Compute -kappa * T 
function sT = getScaledTemperature(rho, u, E, tc)

sT = -tc.gamma * tc.mu / tc.Pr * (E ./ rho - 0.5 * u.^2);

end