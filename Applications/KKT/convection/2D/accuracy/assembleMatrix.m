function [M, S] = assembleMatrix(msh, md, tc, bs, mtol)

% inverse of mass matrix
M = assembleElementMatrix_rect(msh, [0, 0], [0, 0], bs, false, [], mtol);

% element contributions for stiffness matrix
SEx = assembleElementMatrix_rect(msh, [0, 0], [1, 0], bs, false, [], mtol);
SEy = assembleElementMatrix_rect(msh, [0, 0], [0, 1], bs, false, [], mtol);

% face contributions for stiffness matrix
if (any(tc.beta > 0))
    [SFLx, SFLy] = assembleIntFaceMatrix_rect(msh, 'l', bs, md, [], mtol);
end
if (any(tc.beta < 0))
    [SFRx, SFRy] = assembleIntFaceMatrix_rect(msh, 'r', bs, md, [], mtol);
end
if (tc.beta(1) > 0)
    SFx = SFLx;
elseif (tc.beta(1) < 0)
    SFx = SFRx;
end
if (tc.beta(2) > 0)
    SFy = SFLy;
elseif (tc.beta(2) < 0)
    SFy = SFRy;
end

% the global stiffness matrix
S = mychop(tc.beta(1) * (SFx - SEx) + tc.beta(2) * (SFy - SEy), mtol);

end