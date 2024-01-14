function [AU, iUt] = computePartOfTimeDerivative2(msh, md, U, flag, flux, A3, A32, ME, IME, bs)

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
    
    % Initialize the row, column and value arrays 
    row = zeros(bs.nb2, msh.nLElems + 4 * md.nIntLFaces);
    col = row;
    v21 = row;
    v22 = row;
    
    % element contributions 
    e = msh.nLElems;
    row(:, 1 : e) = re + (0 : msh.nLElems - 1) * bs.nb;
    col(:, 1 : e) = ce + (0 : msh.nLElems - 1) * bs.nb;
    
    v21(:, 1 : e) = -reshape(IME * bs.phitw{2} * bs.phi{1} * udr, [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems);
    v22(:, 1 : e) = -reshape(IME * bs.phitw{2} * bs.phi{1} * udm, [bs.nb2, msh.nLElems]) ./ msh.elemJac(:, msh.LElems);
    
    % face contributions  
    s = e + 1; e = e + 4 * nf;
    row(:, s : e) = [re + (leLIDs - 1) * bs.nb, re + (leLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb, re + (reLIDs - 1) * bs.nb];   
    col(:, s : e) = [ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb, ce + (leLIDs - 1) * bs.nb, ce + (reLIDs - 1) * bs.nb];
    
    udrl = bs.phi_face{1, 2} * udr(:, (1 - bs.nb : 0)' + bs.nb * leLIDs);
    udrr = bs.phi_face{1, 1} * udr(:, (1 - bs.nb : 0)' + bs.nb * reLIDs);
    udml = bs.phi_face{1, 2} * udm(:, (1 - bs.nb : 0)' + bs.nb * leLIDs);
    udmr = bs.phi_face{1, 1} * udm(:, (1 - bs.nb : 0)' + bs.nb * reLIDs);

    templ = repelem(ul, 1, bs.nb) .* repmat(bs.phi_face{1, 2}, [1, nf]) + repelem(rhol, 1, bs.nb) .* udrl;
    tempr = repelem(ur, 1, bs.nb) .* repmat(bs.phi_face{1, 1}, [1, nf]) + repelem(rhor, 1, bs.nb) .* udrr;
    v11           = [-(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                     -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                      (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)]; 
    v21(:, s : e) = [ (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * udrl,  [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * udrr,  [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                     -(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * udrl,  [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                     -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * udrr,  [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)]; 
        
    
    templ = repelem(rhol, 1, bs.nb) .* udml;
    tempr = repelem(rhor, 1, bs.nb) .* udmr;
    v12           = [-(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                     -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * templ, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                      (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * tempr, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)]; 
    v22(:, s : e) = [ (0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * udml, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                      (0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 2} * udmr, [bs.nb2, nf]) ./ msh.elemJac(:, leIDs), ...
                     -(0.5 - flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * udml, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs), ...
                     -(0.5 + flux.alpha1) * reshape(IME * bs.phitw_face{1, 1} * udmr, [bs.nb2, nf]) ./ msh.elemJac(:, reIDs)];
               
    % Assemble
    AU = [sparse(row(:, msh.nLElems + 1 : end), col(:, msh.nLElems + 1 : end), v11, NDOF, NDOF), sparse(row(:, msh.nLElems + 1 : end), col(:, msh.nLElems + 1 : end), v12, NDOF, NDOF); A32 * sparse(row, col, v21, NDOF, NDOF), A32 * sparse(row, col, v22, NDOF, NDOF)];
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
    
    % Add contributions of viscous term
    iUt(end / 2 + 1 : end) = A3 * V(:);
else
    iUt = [];
end

end


