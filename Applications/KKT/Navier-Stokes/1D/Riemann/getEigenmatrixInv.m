% the inverse of eigenmatrix 
function RI = getEigenmatrixInv(u, c, tc)

RI = [(tc.gamma - 1) / 4 * u.^2 ./ c.^2 + u ./ (2 * c), -(tc.gamma - 1) / 2 * u ./ c.^2 - 1 / 2 / c, (tc.gamma - 1) / 2 / c.^2; ...
                 1 - (tc.gamma - 1) / 2 * u.^2 ./ c.^2,                  (tc.gamma - 1) * u ./ c.^2,    -(tc.gamma - 1) / c.^2; ...
      (tc.gamma - 1) / 4 * u.^2 ./ c.^2 - u ./ (2 * c), -(tc.gamma - 1) / 2 * u ./ c.^2 + 1 / 2 / c, (tc.gamma - 1) / 2 / c.^2];

end