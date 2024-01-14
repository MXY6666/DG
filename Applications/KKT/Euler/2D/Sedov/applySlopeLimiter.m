% lim = 0 :  no limiter
% lim = 1 : TVD limiter
% lim = 2 : TVB limiter
function U = applySlopeLimiter(msh, md, U, tc, bs, slo)

if (bs.type == 200)
    return
end

switch slo
    
    % no limiter
    case 0
        
        return 
                   
    % slope limiter(TVD or TVB)
    case {1, 2} 
        
        if (mod(floor(bs.type / 10), 10) ~= 0)
            error('Use orthogonal basis functions instead when applying TVD or TVB limiters')
        end
        
        % Set up a small number
        eps = 1e-7; 
              
        U = reshape(U, [bs.nb, 4 * msh.nLElems]);
        U = U(:, (0 : 3)' * msh.nLElems + (1 : msh.nLElems));
        % Compute the difference of mean values in neighboring elements
        um = reshape(U(1, :), [4, msh.nLElems]);
        fdum1 = zeros(4, msh.nLElems);
        bdum1 = zeros(4, msh.nLElems);
        fdum2 = zeros(4, msh.nLElems);
        bdum2 = zeros(4, msh.nLElems);    
        
        % internal faces
        faceIDs = md.intLFaces{1, 3};
        leLIDs  = msh.faceElems(1, faceIDs);
        reLIDs  = msh.faceElems(2, faceIDs);
        fdum1(:, leLIDs) = um(:, reLIDs) - um(:, leLIDs);
        bdum1(:, reLIDs) = fdum1(:, leLIDs);
        
        faceIDs = md.intLFaces{2, 3};
        leLIDs = msh.faceElems(1, faceIDs);
        reLIDs = msh.faceElems(2, faceIDs);
        fdum2(:, leLIDs) = um(:, reLIDs) - um(:, leLIDs);
        bdum2(:, reLIDs) = fdum2(:, leLIDs);  
        
        % boundary faces
        faceIDs = md.bndLFaces{1, 1};
        leLIDs  = msh.faceElems(1, faceIDs);
        nf      = length(faceIDs);
        bdum1(:, leLIDs) = [zeros(1, nf); 2 * um(2, leLIDs); zeros(2, nf)];
        
        faceIDs = md.bndLFaces{3, 1};
        leLIDs  = msh.faceElems(1, faceIDs);  
        nf      = length(faceIDs);
        bdum2(:, leLIDs) = [zeros(2, nf); 2 * um(3, leLIDs); zeros(1, nf)];
        
        % Compute the modified slope of the P1 part, 
        % i.e. the modified second and third coefficients
        ux = reshape(U(2, :), [4, msh.nLElems]);
        uy = reshape(U(3, :), [4, msh.nLElems]);
        ux_mod = zeros(4, msh.nLElems); 
        uy_mod = zeros(4, msh.nLElems);
        if (slo == 1)
            parfor ie = 1 : msh.nLElems  
                % Compute the eigenvector matrix and its inverse 
                temp = um(:, ie);
                rhom = temp(1); m1m = temp(2); m2m = temp(3); Em = temp(4);
        
                pm   = getPressure(rhom, m1m, m2m, Em, tc);
                q2m  = getSquareSpeed(m1m / rhom, m2m / rhom);
                cm   = getSpeedOfSound(rhom, pm, tc);
                Hm   = getEnthalpy(rhom, Em, pm);
                Rmx  = getEigenmatrix_x(m1m / rhom, m2m / rhom, q2m, cm, Hm);
                RImx = getEigenmatrixInv_x(m1m / rhom, m2m / rhom, q2m, cm, tc);
                Rmy  = getEigenmatrix_y(m1m / rhom, m2m / rhom, q2m, cm, Hm);
                RImy = getEigenmatrixInv_y(m1m / rhom, m2m / rhom, q2m, cm, tc); 

                ux_mod(:, ie) = Rmx * minmod(RImx * ux(:, ie), RImx * fdum1(:, ie), RImx * bdum1(:, ie));
                uy_mod(:, ie) = Rmy * minmod(RImy * uy(:, ie), RImy * fdum2(:, ie), RImy * bdum2(:, ie));     
            end
        else
            M = 50;
            hx = msh.elemLength(1, msh.LElems);
            hy = msh.elemLength(2, msh.LElems);
            parfor ie = 1 : msh.nLElems  
                % Compute the eigenvector matrix and its inverse 
                temp = um(:, ie);
                rhom = temp(1); m1m = temp(2); m2m = temp(3); Em = temp(4);
                
                pm   = getPressure(rhom, m1m, m2m, Em, tc);
                q2m  = getSquareSpeed(m1m / rhom, m2m / rhom);
                cm   = getSpeedOfSound(rhom, pm, tc);
                Hm   = getEnthalpy(rhom, Em, pm);
                Rmx  = getEigenmatrix_x(m1m / rhom, m2m / rhom, q2m, cm, Hm);
                RImx = getEigenmatrixInv_x(m1m / rhom, m2m / rhom, q2m, cm, tc);
                Rmy  = getEigenmatrix_y(m1m / rhom, m2m / rhom, q2m, cm, Hm);
                RImy = getEigenmatrixInv_y(m1m / rhom, m2m / rhom, q2m, cm, tc); 

                ux_mod(:, ie) = Rmx * minmod_bar(RImx * ux(:, ie), RImx * fdum1(:, ie), RImx * bdum1(:, ie), M, hx(ie));
                uy_mod(:, ie) = Rmy * minmod_bar(RImy * uy(:, ie), RImy * fdum2(:, ie), RImy * bdum2(:, ie), M, hy(ie));     
            end            
        end
        ux_mod = ux_mod(:)';
        uy_mod = uy_mod(:)';    
            
        % reconstruction 
        isModifiedx = abs(ux_mod - U(2, :)) >= eps;
        isModifiedy = abs(uy_mod - U(3, :)) >= eps;
        U(2, isModifiedx) = ux_mod(isModifiedx);
        U(3, isModifiedy) = uy_mod(isModifiedy);
        if (bs.deg ~= 1)           
            U(4 : end, isModifiedx) = 0;
            U(4 : end, isModifiedy) = 0;            
        end   
        
        U = U(:, (1 : 4 : end)' + (0 : 3));
        U = U(:);
        
        return
end
                                                           
end





