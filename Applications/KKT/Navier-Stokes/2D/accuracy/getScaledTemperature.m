% Compute -kappa * T 
function sT = getScaledTemperature(rho, u, v, E, tc)

sT = -tc.gamma * tc.mu / tc.Pr * (E ./ rho - 0.5 * (u.^2 + v.^2));

end