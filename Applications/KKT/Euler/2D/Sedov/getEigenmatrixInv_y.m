% the inverse of eigenmatrix in y direction 
function RI = getEigenmatrixInv_y(u, v, q2, c, tc)


RI = [(tc.gamma - 1) / 4 * q2 ./ c.^2 + v ./ (2 * c), -(tc.gamma - 1) / 2 * u ./ c.^2, -(tc.gamma - 1) / 2 * v ./ c.^2 - 1 / 2 / c, (tc.gamma - 1) / 2 / c.^2; ...
                 1 - (tc.gamma - 1) / 2 * q2 ./ c.^2,      (tc.gamma - 1) * u ./ c.^2,                  (tc.gamma - 1) * v ./ c.^2,    -(tc.gamma - 1) / c.^2; ...
                                                  -u,                               1,                                           0,                         0; ...
      (tc.gamma - 1) / 4 * q2 ./ c.^2 - v ./ (2 * c), -(tc.gamma - 1) / 2 * u ./ c.^2, -(tc.gamma - 1) / 2 * v ./ c.^2 + 1 / 2 / c, (tc.gamma - 1) / 2 / c.^2];

end