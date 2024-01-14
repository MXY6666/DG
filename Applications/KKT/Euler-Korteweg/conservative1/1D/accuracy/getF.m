function f = getF(rho, u, tc)

f = 0.5 * u.^2 + getFreeEnergyDerivative(rho, tc);

end