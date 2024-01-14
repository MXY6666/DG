% Compute the derivative of free energy density for van der Waals fluids
function dF = getFreeEnergyDerivative(rho, tc)

dF = tc.R * tc.T .* (tc.b ./ (tc.b - rho) + log(rho ./ (tc.b - rho))) - 2 * tc.a * rho;

end