% Compute the DG residual
function res = computeDGResidual(msh, md, U, tc, bs)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(U, [bs.nb, 3 * msh.nLElems]);
U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
u = bs.phi{1} * U;
p = getPressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
if (bs.type == 100)
    res = zeros(1, 3 * msh.nLElems);
else 
    F = getF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    res = -bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************   
%               internal face contributions for DG residual
%**************************************************************************  
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, nf]);
ur = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * reLIDs), [3, nf]);
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = getPressure(ur(1, :), ur(2, :), ur(3, :), tc);
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl); 
Fr = getF(ur(1, :), ur(2, :), ur(3, :), pr);

% dissipation coefficient 
alpha = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
res(:, (-2 : 0)' + 3 * leLIDs) = res(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * reshape(F_hat, [1, 3 * nf]);
res(:, (-2 : 0)' + 3 * reLIDs) = res(:, (-2 : 0)' + 3 * reLIDs) - bs.phitw_face{1, 1} * reshape(F_hat, [1, 3 * nf]);

%**************************************************************************   
%               boundary face contributions for DG residual
%**************************************************************************  
% left boundary
faceIDs = md.bndLFaces{1, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 1} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoL; tc.mL; tc.EL];
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pL;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl);    
Fr = getF(tc.rhoL, tc.mL, tc.EL, tc.pL);

% dissipation coefficient 
alpha = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));

F_hat = 0.5 * (-Fl - Fr - alpha .* (ur - ul));
res(:, (-2 : 0)' + 3 * leLIDs) = res(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 1} * F_hat';

% right boundary
faceIDs = md.bndLFaces{2, 1};
leLIDs  = msh.faceElems(1, faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-2 : 0)' + 3 * leLIDs), [3, 1]);
ur = [tc.rhoR; tc.mR; tc.ER];
pl = getPressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = tc.pR;
Fl = getF(ul(1, :), ul(2, :), ul(3, :), pl);         
Fr = getF(tc.rhoR, tc.mR, tc.ER, tc.pR);

% dissipation coefficient 
alpha = max(getEigenMax(ul(1, :), ul(2, :) ./ ul(1, :), pl, tc), getEigenMax(ur(1, :), ur(2, :) ./ ur(1, :), pr, tc));

F_hat = 0.5 * (Fl + Fr - alpha .* (ur - ul));
res(:, (-2 : 0)' + 3 * leLIDs) = res(:, (-2 : 0)' + 3 * leLIDs) + bs.phitw_face{1, 2} * F_hat';

%**************************************************************************   
%                            DG residual 
%************************************************************************** 
res = res(:, (1 : 3 : end)' + (0 : 2));
res = res(:);

end



