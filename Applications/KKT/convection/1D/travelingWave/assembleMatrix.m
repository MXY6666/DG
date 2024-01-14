function [M, S] = assembleMatrix(msh, md, tc, bs, mtol)

% inverse of mass matrix
M = assembleElementMatrix_line(msh, 0, 0, bs, false, [], mtol);

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% face contributions for stiffness matrix
if (tc.beta > 0)
    SF = assembleIntFaceMatrix_line(msh, 'l', bs, md, [], mtol);
elseif (tc.beta < 0)
    SF = assembleIntFaceMatrix_line(msh, 'r', bs, md, [], mtol);
end

% the global matrix
S = mychop(tc.beta * (SF - SE), mtol);

end