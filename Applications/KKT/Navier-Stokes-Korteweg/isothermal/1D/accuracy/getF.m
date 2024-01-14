function f = getF(rho, u, tc)

f = getFreeEnergyDerivative(rho, tc) - 0.5 * u.^2;

end