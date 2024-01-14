function [IME, M, pr12, pc12, pv12] = assembleMatrix(msh, md, tc, bs, mtol)

% mass matrix at reference line
ME = computeElementMatrix_refLine(0, 0, bs, mtol);

% inverse of mass matrix at reference line
IME = mychop(inv(ME), mtol);

% global mass matrix
[re, ce, ve] = find(ME);
r = re + (0 : 3 * msh.nLElems - 1) * bs.nb;
c = ce + (0 : 3 * msh.nLElems - 1) * bs.nb;
v = repmat(msh.elemJac(:, msh.LElems) .* ve, [1, 3]);
M = sparse(r, c, v, 3 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb); 

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% face contributions for stiffness matrix
SIF = assembleIntFaceMatrix_line(msh, 'c', bs, md, [], mtol);

% boundary face contributions for stiffness matrix
SBF = assembleBndFaceMatrix_line(msh, 'l', bs, md, [], 2, mtol);

% stiffness matrix
S = -SE + SIF + SBF;

% part of block-(1, 2) for the Schur completement
[r, c, v] = find(S);
pr12 = [msh.nLElems * bs.nb + r; 2 * msh.nLElems * bs.nb + r];
pc12 = [c; msh.nLElems * bs.nb + c];
pv12 = [-4 / 3 * tc.mu * v; v];

end

