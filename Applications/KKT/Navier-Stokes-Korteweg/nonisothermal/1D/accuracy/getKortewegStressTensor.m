% Korteweg stress tensor
function xi = getKortewegStressTensor(rho, w, r, tc)

xi = 1 / tc.We * (rho .* r - 0.5 * w.^2);

end