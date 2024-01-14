function alpha = getEigenMax(msh, U, bs)

% Determine the maximum absolute value of f'(uh) = uh on both sides of all faces
U = reshape(U, [bs.nb, msh.nLElems]);
alpha = max(max(abs(bs.phi_face{1} * U), abs(bs.phi_face{2} * U)));

end