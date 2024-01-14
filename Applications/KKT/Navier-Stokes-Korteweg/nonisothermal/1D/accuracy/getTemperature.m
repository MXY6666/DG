function T = getTemperature(rho, u, E, w, tc)

T = 3 / (8 * tc.Cv) * (3 * rho - 0.5 * u.^2 + (E - 0.5 / tc.We * w.^2) ./ rho);

end