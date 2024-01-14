% temperature
function T = getTemperature(rho, u, v, E, w1, w2, tc)

T = 3 / (8 * tc.Cv) * (3 * rho - 0.5 * (u.^2 + v.^2) + (E - 0.5 / tc.We * (w1.^2 + w2.^2)) ./rho);

end