% the eigenmatrix whose columns are the right eigenvectors
function R = getEigenmatrix(u, c, H)

R = [        1,        1,          1; ...
         u - c,        u,      u + c; ...
    H - u .* c, u.^2 / 2, H + u .* c];

end