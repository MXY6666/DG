function Q = computeAuxiliaryVariable(msh, md, U, t, A1, IME, tc, bs)

%**************************************************************************   
%                element and  internal face contributions 
%**************************************************************************
Q = reshape(A1 * U(:), [bs.nb, msh.nLElems]);

%**************************************************************************   
%                     boundary face contributions 
%**************************************************************************
% left boundary
faceID = md.bndLFaces{1, 1};
leID   = msh.faceElems(1, faceID);
leLID  = leID;

u = tc.u(tc.dm(1), t);

Q(:, leLID) = Q(:, leLID) - (IME * (bs.phitw_face{1, 1} * u)) ./ msh.elemJac(:, leID);

% right boundary
faceID = md.bndLFaces{2, 1};
leID   = msh.faceElems(1, faceID);
leLID  = leID;

u = tc.u(tc.dm(2), t);

Q(:, leLID) = Q(:, leLID) + (IME * (bs.phitw_face{1, 2} * u)) ./ msh.elemJac(:, leID);

end