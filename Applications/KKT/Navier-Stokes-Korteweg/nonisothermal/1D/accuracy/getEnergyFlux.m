% energy flux
function jE = getEnergyFlux(rho, w, z, tc)

jE = 1 / tc.We .* rho .* w .* z;

end