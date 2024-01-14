% equation of state
function p = getPressure(rho, m, E, tc)

p = (tc.gamma - 1) * (E - 0.5 * m.^2 ./ rho);

end

