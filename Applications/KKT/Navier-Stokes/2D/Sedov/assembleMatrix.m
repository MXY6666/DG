function [IME, M, pr12, pc12, pv12] = assembleMatrix(msh, md, tc, bs, mtol)

% mass matrix at reference square
ME = computeElementMatrix_refSquare([0, 0], [0, 0], bs);

% inverse of mass matrix at reference square
IME = mychop(inv(ME), mtol);

% global mass matrix
[re, ce, ve] = find(ME);
r = re + (0 : 4 * msh.nLElems - 1) * bs.nb;
c = ce + (0 : 4 * msh.nLElems - 1) * bs.nb;
v = repmat(msh.elemJac(:, msh.LElems) .* ve, [1, 4]);
M = sparse(r, c, v, 4 * msh.nLElems * bs.nb, 4 * msh.nLElems * bs.nb);

% element contributions for stiffness matrix
SEx = assembleElementMatrix_rect(msh, [0, 0], [1, 0], bs, false, [], mtol);
SEy = assembleElementMatrix_rect(msh, [0, 0], [0, 1], bs, false, [], mtol);

% internal face contributions for stiffness matrix
[SIFx, SIFy] = assembleIntFaceMatrix_rect(msh, 'c', bs, md, [], mtol);

% boundary face contributions for stiffness matrix
[SBFx1, SBFy1] = assembleBndFaceMatrix_rect(msh, 'l', bs, md, [], 8, mtol);
[SBFx2, SBFy2] = assembleBndFaceMatrix_rect(msh, 'l', bs, md, [], 7, mtol);

% stiffness matrix
Sx2 = -SEx + SIFx + SBFx2;
Sy2 = -SEy + SIFy + SBFy2;
Sx1 = Sx2 + SBFx1;
Sy1 = Sy2 + SBFy1;

% part of block-(1, 2) for the Schur completement
neb = msh.nLElems * bs.nb;
[rx1, cx1, vx1] = find(Sx1);
[ry1, cy1, vy1] = find(Sy1);
[rx2, cx2, vx2] = find(Sx2);
[ry2, cy2, vy2] = find(Sy2);

pr12 = [neb + [rx1; ry1; ry1; rx1]; 2 * neb + [ry1; rx1; rx1; ry1]; 3 * neb + [rx2; ry2]];
pc12 = [cx1; neb + cy1; 2 * neb + cy1; 3 * neb + cx1; cy1; neb + cx1; 2 * neb + cx1; 3 * neb + cy1; 4 * neb + cx2; 5 * neb + cy2];
pv12 = [-4 / 3 * tc.mu * vx1; -tc.mu * vy1; -tc.mu * vy1; 2 / 3 * tc.mu * vx1; 2 / 3 * tc.mu * vy1; -tc.mu * vx1; -tc.mu * vx1; -4 / 3 * tc.mu * vy1; vx2; vy2];

end

