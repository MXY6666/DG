% Compute the second derivative of free energy density for van der Waals fluids
function d2F = getFreeEnergyDerivative2(rho, tc)

d2F = tc.R * tc.T * tc.b^2 ./ rho ./ (tc.b - rho).^2 - 2 * tc.a;

end