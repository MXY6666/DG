% equation of state
function p = getPressure(rho, T)

p = 8 * rho .* T ./ (3 - rho) - 3 * rho.^2;

end

