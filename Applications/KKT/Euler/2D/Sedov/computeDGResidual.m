% Compute the DG residual
function res = computeDGResidual(msh, md, U, tc, bs)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(U, [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
u = bs.phi{1} * U; 
p = getPressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);
if (bs.type == 200)
    res = zeros(1, 4 * msh.nLElems);
else
    F = getF(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p); 
    G = getG(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);  
    res = -bs.phitw{2} * reshape(msh.elemJxix( :, msh.LElems) .* F, [bs.nep, 4 * msh.nLElems]) ...                  
         -bs.phitw{3} * reshape(msh.elemJetay(:, msh.LElems) .* G, [bs.nep, 4 * msh.nLElems]);  
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
% internal vertical face contributions
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * reLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getF(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pr);

% dissipation coefficient 
alpha = max(getEigenMax(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc), getEigenMax(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :) ./ ur(1 : bs.nfp, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
res(:, (-3 : 0)' + 4 * reLIDs) = res(:, (-3 : 0)' + 4 * reLIDs) - bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% internal horizontal face contributions
faceIDs = md.intLFaces{2, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * reLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
pr = getPressure(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getG(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pr);

% dissipation coefficient 
alpha = max(getEigenMax(ul(1 : bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc), getEigenMax(ur(1 : bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :) ./ ur(1 : bs.nfp, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
res(:, (-3 : 0)' + 4 * reLIDs) = res(:, (-3 : 0)' + 4 * reLIDs) - bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

%**************************************************************************   
%               boundary face contributions for DG residual
%************************************************************************** 
% slip solid wall
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = ul;
ur(bs.nfp + 1 : 2 * bs.nfp, :) = -ul(bs.nfp + 1 : 2 * bs.nfp, :);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getF(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pl);

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc);

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{3, 1};
leLIDs  = msh.faceElems(1, faceIDs); 
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 3} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
ur = ul;
ur(2 * bs.nfp + 1 : 3 * bs.nfp, :) = -ul(2 * bs.nfp + 1 : 3 * bs.nfp, :);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 
Fr = getG(ur(1 : bs.nfp, :), ur(bs.nfp + 1 : 2 * bs.nfp, :), ur(2 * bs.nfp + 1 : 3 * bs.nfp, :), ur(3 * bs.nfp + 1 : end, :), pl);

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc);

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% supersonic outflow
faceIDs = md.bndLFaces{2, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 

F_hat = Fl;
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{4, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 

F_hat = Fl;
res(:, (-3 : 0)' + 4 * leLIDs) = res(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

%**************************************************************************   
%                            DG residual 
%************************************************************************** 
res = res(:, (1 : 4 : end)' + (0 : 3));
res = res(:);

end
