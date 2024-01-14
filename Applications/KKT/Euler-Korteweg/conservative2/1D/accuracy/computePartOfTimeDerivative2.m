function [AU, iUt] = computePartOfTimeDerivative2(msh, md, U, flag, flux, ME, IME, bs)

U = reshape(U, [bs.nb, 2 * msh.nLElems]);
% Compute velocity
rho = bs.phi{1} * U(:, 1 : msh.nLElems);
if (flag == 2)
    V = ME * U(:, msh.nLElems + 1 : end);
    for i = 1 : msh.nLElems
        V(:, i) = (bs.phitw{1} * (rho(:, i) .* bs.phi{1})) \ V(:, i);
    end
else
    V   = U(:, msh.nLElems + 1 : end);
    udr = zeros(bs.nb, msh.nLElems * bs.nb);
    udm = zeros(bs.nb, msh.nLElems * bs.nb);
    
    for i = 1 : msh.nLElems
        fdu = bs.phitw{1} * (rho(:, i) .* bs.phi{1});
        udm_ = fdu \ ME;
        udm(:, (i - 1) * bs.nb + 1 : i * bs.nb) = udm_;
        
        V(:, i) = udm_ * V(:, i);
        fdr = bs.phitw{1} * (bs.phi{1} * V(:, i) .* bs.phi{1});
        udr(:, (i - 1) * bs.nb + 1 : i * bs.nb) = -fdu \ fdr;
    end
end

% face contributions
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

if (flag == 1 || flag == 12)
    % number of degree of freedoms for a single variable
    NDOF = msh.nLElems * bs.nb;

    % row and column indices for an element matrix
    re = repmat((1 : bs.nb)', [bs.nb, 1]);
    ce = repelem((1 : bs.nb)', bs.nb);

    % row, column and value arrays 
    row = [re + (leLIDs - 1) * bs.nb, re + (leLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb];   
    col = [ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb, ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb];

    templ = repelem(ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf]) + repelem(rhol, 1, bs.nb) .* (bs.phi_face{1, 2} * udr(:, (1 - bs.nb : 0)' + bs.nb * leLIDs));
    tempr = repelem(ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]) + repelem(rhor, 1, bs.nb) .* (bs.phi_face{1, 1} * udr(:, (1 - bs.nb : 0)' + bs.nb * reLIDs));
    v1 = [-(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
          -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
           (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
           (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)];   
    
    templ = repelem(rhol, 1, bs.nb) .* (bs.phi_face{1, 2} * udm(:, (1 - bs.nb : 0)' + bs.nb * leLIDs));
    tempr = repelem(rhor, 1, bs.nb) .* (bs.phi_face{1, 1} * udm(:, (1 - bs.nb : 0)' + bs.nb * reLIDs));
    v2 = [-(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
          -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
           (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
           (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)]; 
               
    % Assemble
    AU = [sparse(row, col, v1, NDOF, NDOF), sparse(row, col, v2, NDOF, NDOF); sparse([], [], [], NDOF, 2 * NDOF)];
else
    AU = [];
end
    
if (flag == 2 || flag == 12)
    rhou_hat = (0.5 - flux.alpha1) * (rhol .* ul) + (0.5 + flux.alpha1) * (rhor .* ur);
    iUt = zeros(bs.nb, 2 * msh.nLElems);
    iUt(:, leLIDs) = iUt(:, leLIDs) - bs.phitw_face{1, 2} * rhou_hat;
    iUt(:, reLIDs) = iUt(:, reLIDs) + bs.phitw_face{1, 1} * rhou_hat;

    % Take care of the mass matrix
    iUt(:, 1 : msh.nLElems) = (IME * iUt(:, 1 : msh.nLElems)) ./ msh.elemJac(:, msh.LElems);
    iUt = iUt(:);
    
else
    iUt = [];
end

end


