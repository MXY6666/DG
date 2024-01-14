function [AU, iUt] = computePartOfTimeDerivative1(msh, md, U, t, flux, AE, A2, ME, IME, tc, quad1, bs)

U = reshape(U, [bs.nb, 2 * msh.nLElems]);
% number of degree of freedoms for a single variable
ndof = msh.nLElems * bs.nb;

% row and column indices for an element matrix
re = repmat((1 : bs.nb)', [bs.nb, 1]);
ce = repelem((1 : bs.nb)', bs.nb);

% Initialize the row, column and value arrays for part of the bottom left 
% matrix 
r21 = zeros(bs.nb2, msh.nLElems + 4 * md.nIntLFaces);
c21 = r21;
v21 = r21;

% Evaluate some quantities
rho = bs.phi{1} * U(:, 1 : msh.nLElems);
%**************************************************************************
%                   Compute some auxiliary variables
%**************************************************************************
% Compute velocity
V = ME * U(:, msh.nLElems + 1 : end);
for i = 1 : msh.nLElems
    V(:, i) = (bs.phitw{1} * (rho(:, i) .* bs.phi{1})) \ V(:, i);
end
u = bs.phi{1} * V;

% Compute F 
F = IME * (bs.phitw{1} * getF(rho, u, tc));

%**************************************************************************
%                         element contributions   
%**************************************************************************
iUt(1 : bs.nb, msh.nLElems + 1 : 2 * msh.nLElems) = bs.phitw{2} * (rho .* u.^2) - bs.phitw{1} * (rho .* (bs.phi{2} * F + u .* (bs.phi{2} * V)));

e = msh.nLElems;
r21(:, 1 : e) = re + (0 : msh.nLElems - 1) * bs.nb;
c21(:, 1 : e) = ce + (0 : msh.nLElems - 1) * bs.nb;
v21(:, 1 : e) = reshape(IME * bs.phitw{1} * (repelem(rho, 1, bs.nb) .* repmat(bs.phi{2}, [1, msh.nLElems])), [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems);

%**************************************************************************
%                         face contributions   
%**************************************************************************
faceIDs = md.intLFaces{1, 3};
leIDs   = msh.faceElems(1, faceIDs);
reIDs   = msh.faceElems(2, faceIDs);
leLIDs  = leIDs;
reLIDs  = reIDs;
nf      = length(faceIDs);

rhol = bs.phi_face{1, 2} * U(:, leLIDs); 
rhor = bs.phi_face{1, 1} * U(:, reLIDs); 
ul   = bs.phi_face{1, 2} * V(:, leLIDs); 
ur   = bs.phi_face{1, 1} * V(:, reLIDs);
fl   = bs.phi_face{1, 2} * F(:, leLIDs); 
fr   = bs.phi_face{1, 1} * F(:, reLIDs);

rhou_hat  = (0.5 - flux.alpha1) * (rhol .* ul) + (0.5 + flux.alpha1) * (rhor .* ur);
rhou2_hat = 0.5 * (ul + ur) .* rhou_hat;
fu2_jump = (fr + 0.5 * ur.^2) - (fl + 0.5 * ul.^2);
iUt(:, msh.nLElems + leLIDs) = iUt(:, msh.nLElems + leLIDs) - bs.phitw_face{1, 2} * rhou2_hat - (0.5 - flux.alpha1) * bs.phitw_face{1, 2} * (rhol .* fu2_jump);
iUt(:, msh.nLElems + reLIDs) = iUt(:, msh.nLElems + reLIDs) + bs.phitw_face{1, 1} * rhou2_hat - (0.5 + flux.alpha1) * bs.phitw_face{1, 1} * (rhor .* fu2_jump);

s = e + 1; e = e + 4 * nf;
r21(:, s : e) = [re + (leLIDs - 1) * bs.nb, re + (leLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb];
c21(:, s : e) = [ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb, ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb];              
v21(:, s : e) = [-(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * (repelem(rhol, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                  (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * (repelem(rhol, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                 -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * (repelem(rhor, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                  (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * (repelem(rhor, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf])), [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)];

%**************************************************************************   
%                     source term contributions 
%**************************************************************************
[src1, src2] = getSourceTerm(msh.elemCenter(:, msh.LElems) + 0.5 * msh.elemLength(:, msh.LElems) .* quad1.points, t, tc);
iUt(:, 1 : msh.nLElems)       = iUt(:, 1 : msh.nLElems)       + msh.elemJac(:, msh.LElems) .* (bs.phitw{1} * src1);
iUt(:, msh.nLElems + 1 : end) = iUt(:, msh.nLElems + 1 : end) + msh.elemJac(:, msh.LElems) .* (bs.phitw{1} * src2);

%************************************************************************** 
%                            post-processing
%**************************************************************************
% Take care of the mass matrix
iUt = (IME * iUt) ./ repmat(msh.elemJac(:, msh.LElems), [1, 2]);
iUt = iUt(:);

%**************************************************************************
%                         global matrix   
%**************************************************************************
AU = [sparse([], [], [], ndof, ndof), -AE; sparse(r21, c21, v21, ndof, ndof) * A2, sparse([], [], [], ndof, ndof)];

end



