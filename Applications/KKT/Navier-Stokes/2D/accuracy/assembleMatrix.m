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

% face contributions for stiffness matrix
[SFx, SFy] = assembleIntFaceMatrix_rect(msh, 'c', bs, md, [], mtol);

% stiffness matrix
Sx = -SEx + SFx;
Sy = -SEy + SFy;

% part of block-(1, 2) for the Schur completement
neb = msh.nLElems * bs.nb;
[rx, cx, vx] = find(Sx);
[ry, cy, vy] = find(Sy);

pr12 = [neb + [rx; ry; ry; rx]; 2 * neb + [ry; rx; rx; ry]; 3 * neb + [rx; ry]];
pc12 = [cx; neb + cy; 2 * neb + cy; 3 * neb + cx; cy; neb + cx; 2 * neb + cx; 3 * neb + cy; 4 * neb + cx; 5 * neb + cy];
pv12 = [-4 / 3 * tc.mu * vx; -tc.mu * vy; -tc.mu * vy; 2 / 3 * tc.mu * vx; 2 / 3 * tc.mu * vy; -tc.mu * vx; -tc.mu * vx; -4 / 3 * tc.mu * vy; vx; vy];

end

