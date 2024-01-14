function [AU, iUt] = computeTimeDerivative(msh, U, t, A11, A12, Acap, IME, tc, quad1, bs)

% number of degree of freedoms for a single variable
ndof = msh.nLElems * bs.nb;

% row and column indices for an element matrix
re = repmat((1 : bs.nb)', [bs.nb, 1]);
ce = repelem((1 : bs.nb)', bs.nb);

% Evaluate values 
u = bs.phi{1} * reshape(U, [bs.nb, 2 * msh.nLElems]);
%**************************************************************************   
%         Compute some matrices associated with auxiliary variables
%**************************************************************************
row = re + (0 : msh.nLElems - 1) * bs.nb;
col = ce + (0 : msh.nLElems - 1) * bs.nb;
vm  = reshape(IME * bs.phitw{1} * (repelem(u(:, 1 : msh.nLElems), 1, bs.nb) .* repmat(bs.phi{1}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]);

%**************************************************************************
%                    first order term contributions   
%**************************************************************************
F = IME * (bs.phitw{1} * getF(u(:, 1 : msh.nLElems), u(:, msh.nLElems + 1 : end), tc));
iUt(ndof + 1 : 2 * ndof, 1) = -(A11 * F(:));

%**************************************************************************   
%                       source term contributions 
%**************************************************************************
[src1, src2] = getSourceTerm(msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points, t, tc);
iUt(1 : ndof)       =                       reshape(IME * (bs.phitw{1} * src1), [ndof, 1]);
iUt(ndof + 1 : end) = iUt(ndof + 1 : end) + reshape(IME * (bs.phitw{1} * src2), [ndof, 1]);

%**************************************************************************
%                         global matrix   
%**************************************************************************
em = sparse([], [], [], ndof, ndof);
AU = [em, -A12 * sparse(row, col, vm, ndof, ndof); Acap, em];

end

