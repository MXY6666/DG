% the flux function in y direction
function G = getG(rho, m1, m2, E, p)

G = [m2; m1 .* m2 ./ rho; m2.^2 ./ rho + p; m2 ./ rho .* (E + p)];

end