% lim = 0 :  no limiter
% lim = 1 : TVD limiter
% lim = 2 : TVB limiter
function U = applySlopeLimiter(msh, md, U, tc, bs, slo)

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
                        
        U = reshape(U, [bs.nb, 3 * msh.nLElems]);
        U = U(:, (0 : 2)' * msh.nLElems + (1 : msh.nLElems));
        % Compute the difference of mean values in neighboring elements
        um   = reshape(U(1, :), [3, msh.nLElems]);
        fdum = zeros(3, msh.nLElems);
        bdum = zeros(3, msh.nLElems);
        
        faceIDs = md.intLFaces{1, 3};
        leLIDs  = msh.faceElems(1, faceIDs);
        reLIDs  = msh.faceElems(2, faceIDs);
        fdum(:, leLIDs) = um(:, reLIDs) - um(:, leLIDs);
        bdum(:, reLIDs) = fdum(:, leLIDs); 
        
        faceID = md.bndLFaces{1, 1};  
        leLID = msh.faceElems(1, faceID);
        bdum(:, leLID) = um(:, leLID) - [tc.rhoL; tc.mL; tc.EL]; 

        faceID = md.bndLFaces{2, 1};
        leLID = msh.faceElems(1, faceID);
        fdum(:, leLID) = [tc.rhoR; tc.mR; tc.ER] - um(:, leLID); 
                 
        % Compute the modified slop of P1 part, 
        % i.e. the modified second coefficient
        ux     = reshape(U(2, :), [3, msh.nLElems]);
        ux_mod = zeros(3, msh.nLElems);
        if (slo == 1)
            for ie = 1 : msh.nLElems        
                % Compute the eigenmatrix and its inverse 
                pm  = getPressure(um(1, ie), um(2, ie), um(3, ie), tc);
                cm  = getSpeedOfSound(um(1, ie), pm, tc);
                Hm  = getEnthalpy(um(1, ie), um(3, ie), pm);
                Rm  = getEigenmatrix(um(2, ie) / um(1, ie), cm, Hm);
                RIm = getEigenmatrixInv(um(2, ie) / um(1, ie), cm, tc);
                
                ux_mod(:, ie) = Rm * minmod(RIm * ux(:, ie), RIm * fdum(:, ie), RIm * bdum(:, ie));
            end  
        else
            if (tc.pro == 9)
                M = [10; 1; 10];
            else
                M = [50; 50; 50];
            end
            h = msh.elemLength(:, msh.LElems);
            for ie = 1 : msh.nLElems        
                % Compute the eigenmatrix and its inverse 
                pm  = getPressure(um(1, ie), um(2, ie), um(3, ie), tc);
                cm  = getSpeedOfSound(um(1, ie), pm, tc);
                Hm  = getEnthalpy(um(1, ie), um(3, ie), pm);
                Rm  = getEigenmatrix(um(2, ie) / um(1, ie), cm, Hm);
                RIm = getEigenmatrixInv(um(2, ie) / um(1, ie), cm, tc);
                
                ux_mod(:, ie) = Rm * minmod_bar(RIm * ux(:, ie), RIm * fdum(:, ie), RIm * bdum(:, ie), M, h(1, ie));
            end                       
        end
        ux_mod = ux_mod(:)';
        
        % Reconstruction
        isModified = abs(ux_mod - U(2, :)) >= eps;
        U(2, isModified) = ux_mod(isModified);
        if (bs.deg ~= 1)           
            U(3 : end, isModified) = 0;            
        end
       
        U = U(:, (1 : 3 : end)' + (0 : 2));
        U = U(:);
        
        return
end
                                                           
end

