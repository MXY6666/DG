function [IME, A, M, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i] = assembleMatrix(msh, md, tc, bs, mtol)

% mass matrix at reference line
ME = computeElementMatrix_refLine(0, 0, bs, mtol);

% inverse of mass matrix at reference line
IME = mychop(inv(ME), mtol);

% inverse of mass matrix for a singal variable
IM = assembleElementMatrix_line(msh, 0, 0, bs, true, IME, mtol);

% global mass matrix
[re, ce, ve] = find(ME);
r = re + (0 : 3 * msh.nLElems - 1) * bs.nb;
c = ce + (0 : 3 * msh.nLElems - 1) * bs.nb;
v = repmat(msh.elemJac(:, msh.LElems) .* ve, [1, 3]);
M = sparse(r, c, v, 3 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb); 

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% face contributions for stiffness matrix
SFC = assembleIntFaceMatrix_line(msh, 'c', bs, md, [], mtol);

% stiffness matrix
S = -SE + SFC;        

% Multiply inverse of mass matrix and stiffness matrix
A = mychop(IM * S, mtol);       

% part of block-(1, 2) for the Schur completement  
[r, c, v] = find(S);
pr12 = [msh.nLElems * bs.nb + r; 2 * msh.nLElems * bs.nb + r];
pc12 = [2 * msh.nLElems * bs.nb + c; 3 * msh.nLElems * bs.nb + c];
pv12 = [-4 / 3 / tc.Re * v; v];  

% part of block-(2, 1) for the Schur completement 
[r, c, v] = find(A);
pr21 = r;
pc21 = c;
pv21 = v;

% part of inverse of block-(2, 2) for the Schur completement 
pr22i = [(1 : 4 * msh.nLElems * bs.nb)'; msh.nLElems * bs.nb + r];
pc22i = [(1 : 4 * msh.nLElems * bs.nb)'; c];
pv22i = [-ones(4 * msh.nLElems * bs.nb, 1); -v];
        
end

