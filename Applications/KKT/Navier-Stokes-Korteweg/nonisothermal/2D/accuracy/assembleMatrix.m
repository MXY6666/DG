function [IME, Ax, Ay, M, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i] = assembleMatrix(msh, md, tc, bs, mtol)

% mass matrix at reference square
ME = computeElementMatrix_refSquare([0, 0], [0, 0], bs);

% inverse of mass matrix at reference square
IME = mychop(inv(ME), mtol);

% inverse of mass matrix for a singal variable
IM = assembleElementMatrix_rect(msh, [0, 0], [0, 0], bs, true, IME, mtol);

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

% Multiply inverse of mass matrix and stiffness matrix
Ax = mychop(IM * Sx, mtol);
Ay = mychop(IM * Sy, mtol); 
            
neb = msh.nLElems * bs.nb;
% part of block-(1, 2) for the Schur completement
[rx, cx, vx] = find(Sx);
[ry, cy, vy] = find(Sy);

pr12 = [neb + [rx; ry; ry; rx]; 2 * neb + [ry; rx; rx; ry]; 3 * neb + [rx; ry]];
pc12 = [3 * neb + cx; 4 * neb + cy; 5 * neb + cy; 6 * neb + cx; 3 * neb + cy; 4 * neb + cx; 5 * neb + cx; 6 * neb + cy; 7 * neb + cx; 8 * neb + cy];
pv12 = [-4 / 3 / tc.Re * vx; -1 / tc.Re * vy; -1 / tc.Re * vy; 2 / 3 / tc.Re * vx; 2 / 3 / tc.Re * vy; -1 / tc.Re * vx; -1 / tc.Re * vx; -4 / 3 / tc.Re * vy; vx; vy];

% part of block-(2, 1) for the Schur completement 
[rx, cx, vx] = find(Ax);
[ry, cy, vy] = find(Ay);      

pr21 = [rx; neb + ry];
pc21 = [cx; cy];
pv21 = [vx; vy];

% part of inverse of block-(2, 2) for the Schur completement 
pr22i = [(1 : 9 * neb)'; 2 * neb + rx; 2 * neb + ry];
pc22i = [(1 : 9 * neb)'; cx; neb + cy];
pv22i = [-ones(9 * neb, 1); -vx; -vy];
    
end
