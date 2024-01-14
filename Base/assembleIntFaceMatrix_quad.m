% msh  : linear quadrilateral mesh in 2D
% flux : numerical flux, 'l' means upwind, 'r'means downwind, 'c' means
%        central
% bs   : basis function set data (a struct) or just the degree of 
%        polynomial (an integer)
% md   : mesh data
% FM   : face matrix at reference geometry
% wind : 2D vector representing wind direction
% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function [Ax, Ay] = assembleIntFaceMatrix_quad(msh, flux, bs, md, FM, wind, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 211) && (msh.type ~= 212)
    error('Wrong mesh type')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_quad(msh);
end 

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference square
if (nargin < 5) || isempty(FM)
    if (msh.type == 211)
        FM = cell(4, 4);
    else
        FM = cell(12, 12);
    end
    for lfn = 1 : 4
        for rfn = 1 : 4
            temp1 = computeFaceMatrix_refSquare(lfn, rfn, 2, 'l', bs, mtol);
            temp2 = computeFaceMatrix_refSquare(lfn, rfn, 2, 'r', bs, mtol);   
            FM{lfn, rfn} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
        end
    end        
    if (msh.type == 212)
        for lfn = 1 : 4
            for rfn = 5 : 12
                temp1 = computeFaceMatrix_refSquare(lfn, rfn, 2, 'l', bs, mtol);
                temp2 = computeFaceMatrix_refSquare(lfn, rfn, 2, 'r', bs, mtol);
                FM{lfn, rfn} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
                
                temp1 = computeFaceMatrix_refSquare(rfn, lfn, 2, 'l', bs, mtol);
                temp2 = computeFaceMatrix_refSquare(rfn, lfn, 2, 'r', bs, mtol);
                FM{rfn, lfn} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
            end
        end            
    end
end

if (nargin < 6) || isempty(wind)
    wind = [1; 1];
end
if (length(wind) ~= 2)
    error('Wrong size of wind')
end

nb = size(FM{1, 1}{1, 1}, 1);
nb2 = nb * nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 0;
re = repmat((1 : nb)', [nb, 1]);
ce = repelem((1 : nb)', nb);
if (msh.type == 211)
    switch flux
        case 'l'
            faceIDs = [md.intLFaces{:, 3}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                J         = msh.faceJac(:, ID);
                
                leLID = faceElems(1); lfn = faceNums(1);
                reLID = faceElems(2); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2);               
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, leLID - 1] * nb;
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [reLID - 1, reLID - 1] * nb;  
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];                    
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e)  = re + [leLID - 1, leLID - 1, reLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, reLID - 1, leLID - 1, reLID - 1] * nb;
                    vx(:, s : e) = 0.5 * normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = 0.5 * normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];                    
                end
            end 
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
        case 'r'
            faceIDs = [md.intLFaces{:, 3}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                J         = msh.faceJac(:, ID);

                leLID = faceElems(1); lfn = faceNums(1);
                reLID = faceElems(2); rfn = faceNums(2);                
                ip    = normalx * wind(1) + normaly * wind(2); 
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [reLID - 1, reLID - 1] * nb;  
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)]; 
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, leLID - 1] * nb;
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];                   
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e)  = re + [leLID - 1, leLID - 1, reLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, reLID - 1, leLID - 1, reLID - 1] * nb;
                    vx(:, s : e) = 0.5 * normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = 0.5 * normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];                    
                end
            end     
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);            
        case 'c'
            r = zeros(4 * md.nIntLFaces * nb2, 1); c = r; vx = r; vy = r;   
            for i = 1 : 16
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.faceElems(1, faceIDs);
                    reLIDs = msh.faceElems(2, faceIDs);
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2});
                    vll = 0.5 * vll;
                    vlr = 0.5 * vlr;
                    vrl = 0.5 * vrl; 
                    vrr = 0.5 * vrr;
                                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e)  = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e)  = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    vx(s : e) = reshape(msh.faceNormalx(:, faceIDs) .* msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr], [l, 1]); 
                    vy(s : e) = reshape(msh.faceNormaly(:, faceIDs) .* msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr], [l, 1]); 
                end              
            end   
            Ax = sparse(r(1 : e), c(1 : e), vx(1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(1 : e), c(1 : e), vy(1 : e), msh.nLElems * nb, msh.nLElems * nb);            
        otherwise
            error('Wrong flux')
    end
else
    switch flux
        case 'l'
            faceIDs = [md.intLFaces{:, 3}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                J         = msh.faceJac(:, ID);                
                
                leLID = msh.elemLID(faceElems(1)); lfn = faceNums(1);
                reLID = msh.elemLID(faceElems(2)); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2);                 
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, leLID - 1] * nb;
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [reLID - 1, reLID - 1] * nb; 
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];                    
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e)  = re + [leLID - 1, leLID - 1, reLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, reLID - 1, leLID - 1, reLID - 1] * nb;
                    vx(:, s : e) = 0.5 * normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = 0.5 * normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];                    
                end
            end
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);            
        case 'r'
            faceIDs = [md.intLFaces{:, 3}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                J         = msh.faceJac(:, ID);                
                
                leLID = msh.elemLID(faceElems(1)); lfn = faceNums(1);
                reLID = msh.elemLID(faceElems(2)); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2);
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [reLID - 1, reLID - 1] * nb; 
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 2}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e)  = re + [leLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, leLID - 1] * nb;
                    vx(:, s : e) = normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];
                    vy(:, s : e) = normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{2, 1}(:)];                   
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e)  = re + [leLID - 1, leLID - 1, reLID - 1, reLID - 1] * nb;
                    c(:, s : e)  = ce + [leLID - 1, reLID - 1, leLID - 1, reLID - 1] * nb;
                    vx(:, s : e) = 0.5 * normalx * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];
                    vy(:, s : e) = 0.5 * normaly * J * [FM{lfn, rfn}{1, 1}(:), FM{lfn, rfn}{1, 2}(:), FM{lfn, rfn}{2, 1}(:), FM{lfn, rfn}{2, 2}(:)];                    
                end
            end  
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);            
        case 'c'  
            r = zeros(4 * md.nIntLFaces * nb2, 1); c = r; vx = r; vy = r;   
            for i = 1 : 80
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    J      = msh.faceJac(:, faceIDs); 
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2});
                    vll = 0.5 * vll;
                    vlr = 0.5 * vlr;
                    vrl = 0.5 * vrl; 
                    vrr = 0.5 * vrr;                    
                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e)  = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e)  = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    vx(s : e) = reshape(msh.faceNormalx(:, faceIDs) .* J .* [vll; vlr; vrl; vrr], [l, 1]); 
                    vy(s : e) = reshape(msh.faceNormaly(:, faceIDs) .* J .* [vll; vlr; vrl; vrr], [l, 1]); 
                end              
            end   
            Ax = sparse(r(1 : e), c(1 : e), vx(1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(1 : e), c(1 : e), vy(1 : e), msh.nLElems * nb, msh.nLElems * nb); 
        otherwise
            error('Wrong flux')
    end
end


end
