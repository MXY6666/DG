% total specific enthalpy
function H = getEnthalpy(rho, E, p)

H = (E + p) ./ rho;

end