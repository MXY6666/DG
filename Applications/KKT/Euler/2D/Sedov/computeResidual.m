% Compute the residual when applying DIRK time integration method
function [res, LStage, waveVel] = computeResidual(msh, md, U0, KStage, Lamb, LStage, dt, stage, M, tc, bs, bt, KKT)

% Compute DG coefficients for current Runge-Kutta stage
U = reshape(getStageSolution(U0, KStage, dt, stage, bt), [bs.nb, 4 * msh.nLElems]);
U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));

%**************************************************************************
%                 element contributions for DG residual
%************************************************************************** 
u = bs.phi{1} * U; 
p = getPressure(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), tc);
if (bs.type == 200)
    Lh = zeros(1, 4 * msh.nLElems);
else
    F = getF(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p); 
    G = getG(u(:, 1 : 4 : end), u(:, 2 : 4 : end), u(:, 3 : 4 : end), u(:, 4 : 4 : end), p);  
    Lh = -bs.phitw{2} * reshape(msh.elemJxix( :, msh.LElems) .* F, [bs.nep, 4 * msh.nLElems]) ...                  
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
waveVel = max(alpha, [], 'all');

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur - ul));
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
Lh(:, (-3 : 0)' + 4 * reLIDs) = Lh(:, (-3 : 0)' + 4 * reLIDs) - bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

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
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = 0.5 * (Fl + Fr - repmat(alpha, [4, 1]) .* (ur - ul));
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);
Lh(:, (-3 : 0)' + 4 * reLIDs) = Lh(:, (-3 : 0)' + 4 * reLIDs) - bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

%**************************************************************************   
%               boundary face contributions for DG residual
%************************************************************************** 
% slip solid wall boundary face contributions
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
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 1} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

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
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = 0.5 * (-Fl - Fr - repmat(alpha, [4, 1]) .* (ur - ul));
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 3} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

% supersonic outflow boundary face contributions
faceIDs = md.bndLFaces{2, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 2} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getF(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc);
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = Fl;
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 2} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

faceIDs = md.bndLFaces{4, 2};
leLIDs  = msh.faceElems(1, faceIDs);
J       = msh.faceJac(:, faceIDs);
nf      = length(faceIDs);

ul = reshape(bs.phi_face{1, 4} * U(:, (-3 : 0)' + 4 * leLIDs), [4 * bs.nfp, nf]);
pl = getPressure(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), tc);
Fl = getG(ul(1 : bs.nfp, :), ul(bs.nfp + 1 : 2 * bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :), ul(3 * bs.nfp + 1 : end, :), pl); 

% dissipation coefficient 
alpha = getEigenMax(ul(1 : bs.nfp, :), ul(2 * bs.nfp + 1 : 3 * bs.nfp, :) ./ ul(1 : bs.nfp, :), pl, tc);
waveVel = max(waveVel, max(alpha, [], 'all'));

F_hat = Fl;
Lh(:, (-3 : 0)' + 4 * leLIDs) = Lh(:, (-3 : 0)' + 4 * leLIDs) + bs.phitw_face{1, 4} * reshape(J .* F_hat, [bs.nfp, 4 * nf]);

%**************************************************************************   
%                          residual vector
%************************************************************************** 
% Compute DG residual for current Runge-Kutta stage
LStage(:, stage.Ls) = reshape(Lh(:, (1 : 4 : end)' + (0 : 3)), [4 * msh.nLElems * bs.nb, 1]);

% DG discretization contribution for residual vector
res = zeros(4 * KKT.NDOF + 2 * KKT.LDOF, 1);
if (bt.form == 1) 
    res(1 : 4 * KKT.NDOF) = M * (KStage(:, 1) - U0) + dt * (LStage(:, 1 : stage.s) * bt.A(stage.s, 1 : stage.s)');
elseif (bt.form == 2)  
    res(1 : 4 * KKT.NDOF) = M * KStage(:, stage.s) + LStage(:, 1);
end

% Add positivity constraints contribution
if (KKT.positivity)
    if (bt.form == 1) 
        cdt = 1;
    elseif (bt.form == 2)
        cdt = dt * bt.diagCoe;       
    end    
    Lamb = reshape(Lamb, [bs.nep, 2 * msh.nLElems]);
    v = [u(:, 2 : 4 : end) ./ u(:, 1 : 4 : end), u(:, 3 : 4 : end) ./ u(:, 1 : 4 : end)];
        
    res(1 : KKT.NDOF) = res(1 : KKT.NDOF) - cdt * reshape(bs.phi{1}' * Lamb(:, 1 : msh.nLElems), [KKT.NDOF, 1]);
    for i = 1 : 4
        dp = getPressureDeriv(v(:, 1 : end / 2), v(:, end / 2 + 1 : end), i, tc);
        res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) = res((i - 1) * KKT.NDOF + 1 : i * KKT.NDOF) - cdt * reshape(bs.phi{1}' * (dp .* Lamb(:, msh.nLElems + 1 : end)), [KKT.NDOF, 1]);
    end

    g = KKT.Umin - [u(:, 1 : 4 : end), p];
    res(4 * KKT.NDOF + 1 : end) = min(-g(:), Lamb(:));
end

end
