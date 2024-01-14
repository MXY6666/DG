% the flux function 
function F = getF(rho, m, E, p)

F = [m; m.^2 ./ rho + p; m ./ rho .* (E + p)];

end

