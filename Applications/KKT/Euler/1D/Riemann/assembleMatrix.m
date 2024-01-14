function M = assembleMatrix(msh, bs, mtol)

% mass matrix at reference line
ME = computeElementMatrix_refLine(0, 0, bs, mtol);

% global mass matrix
[re, ce, ve] = find(ME);
r = re + (0 : 3 * msh.nLElems - 1) * bs.nb;
c = ce + (0 : 3 * msh.nLElems - 1) * bs.nb;
v = repmat(msh.elemJac(:, msh.LElems) .* ve, [1, 3]);
M = sparse(r, c, v, 3 * msh.nLElems * bs.nb, 3 * msh.nLElems * bs.nb);

end