% equation of state
function p = getPressure(rho, m1, m2, E, tc)

p = (tc.gamma - 1) * (E - 1 / 2 * (m1.^2 + m2.^2) ./ rho);

end

