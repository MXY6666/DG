% Compute the free energy density for van der Waals fluids
function F = getFreeEnergy(rho, tc)

F = tc.R * tc.T .* rho .* log(rho ./ (tc.b - rho)) - tc.a * rho.^2;

end