% msh  : linear hexahedral mesh in 3D
% flux : numerical flux, 'l' means upwind, 'r' means downwind, 'c' means
%        central
% bs   : basis function set data (a struct) or just the degree of 
%        polynomial (an integer)
% md   : mesh data
% FM   : face matrix at reference geometry
% wind : 3D vector representing wind direction
% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function [Ax, Ay, Az] = assembleIntFaceMatrix_hexa(msh, flux, bs, md, FM, wind, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 311) && (msh.type ~= 312)
    error('Wrong mesh type')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_hexa(msh);
end 

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference square
% Note, some combinations of f2f, lfn and rfn may not make any sense
if (nargin < 5) || isempty(FM)
    if (msh.type == 311)
        FM = cell(6, 6, 8);
    else
        FM = cell(30, 30, 8);
    end
    for f2f = 1 : 8
        for rfn = 1 : 6
            for lfn = 1 : 6
                temp1 = computeFaceMatrix_refCube(lfn, rfn, f2f, 'l', bs, mtol);
                temp2 = computeFaceMatrix_refCube(lfn, rfn, f2f, 'r', bs, mtol);  
                FM{lfn, rfn, f2f} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
            end
        end      
    end
    if (msh.type == 312)
        for f2f = 1 : 8    
            for rfn = 7 : 30
                for lfn = 1 : 6
                    temp1 = computeFaceMatrix_refCube(lfn, rfn, f2f, 'l', bs, mtol);
                    temp2 = computeFaceMatrix_refCube(lfn, rfn, f2f, 'r', bs, mtol);
                    FM{lfn, rfn, f2f} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
                    
                    temp1 = computeFaceMatrix_refCube(rfn, lfn, f2f, 'l', bs, mtol);
                    temp2 = computeFaceMatrix_refCube(rfn, lfn, f2f, 'r', bs, mtol);
                    FM{rfn, lfn, f2f} = {temp1{1, 1}, temp2{1, 2}; temp1{2, 1}, temp2{2, 2}};
                end
            end            
        end
    end
end

if (nargin < 6) || isempty(wind)
    wind = [1; 1; 1];
end
if (length(wind) ~= 3)
    error('Wrong size of wind')
end

nb = size(FM{1, 1}{1, 1}, 1);
nb2 = nb * nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 0;
re = repmat((1 : nb)', [nb, 1]);
ce = repelem((1 : nb)', nb);
if (msh.type == 311)
    switch flux
        case 'l'
            faceIDs = [md.intLFaces{:, 3, :}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r; vz = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                f2f       = msh.faceR2L(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                normalz   = msh.faceNormalz(:, ID);
                J         = msh.faceJac(:, ID);
                
                leLID = faceElems(1); lfn = faceNums(1);
                reLID = faceElems(2); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2) + normalz * wind(3);  
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (leLID - 1) * nb];
                    v = J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{2, 1}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (reLID - 1) * nb, ce + (reLID - 1) * nb];  
                    v = J * [FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 2}(:)];                 
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (leLID - 1) * nb, re + (reLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (reLID - 1) * nb, ce + (leLID - 1) * nb, ce + (reLID - 1) * nb];
                    v = 0.5 * J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 1}(:), FM{lfn, rfn, f2f}{2, 2}(:)];                 
                end
                vx(:, s : e) = normalx * v;
                vy(:, s : e) = normaly * v;
                vz(:, s : e) = normalz * v;                
            end             
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Az = sparse(r(:, 1 : e), c(:, 1 : e), vz(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
        case 'r'
            faceIDs = [md.intLFaces{:, 3, :}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r; vz = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                f2f       = msh.faceR2L(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                normalz   = msh.faceNormalz(:, ID);
                J         = msh.faceJac(:, ID);
                
                leLID = faceElems(1); lfn = faceNums(1);
                reLID = faceElems(2); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2) + normalz * wind(3); 
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (reLID - 1) * nb, ce + (reLID - 1) * nb];  
                    v = J * [FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 2}(:)]; 
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (leLID - 1) * nb];
                    v = J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{2, 1}(:)];                 
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (leLID - 1) * nb, re + (reLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (reLID - 1) * nb, ce + (leLID - 1) * nb, ce + (reLID - 1) * nb];
                    v = 0.5 * J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 1}(:), FM{lfn, rfn, f2f}{2, 2}(:)];
                end
                vx(:, s : e) = normalx * v;
                vy(:, s : e) = normaly * v;
                vz(:, s : e) = normalz * v;
            end     
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb); 
            Az = sparse(r(:, 1 : e), c(:, 1 : e), vz(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
        case 'c'
            r = zeros(4 * md.nIntLFaces * nb2, 1); c = r; vx = r; vy = r; vz = r;
            i = 1;
            while (i <= 288)
                f2f = ceil(i / 36);
                if (md.nifs(f2f) == 0)
                    i = i + 36;
                    continue;
                end
                j       = i - 36 * (f2f - 1);
                lfn     = md.intLFaces{j, 1, f2f};
                rfn     = md.intLFaces{j, 2, f2f};
                faceIDs = md.intLFaces{j, 3, f2f};
                if ~isempty(faceIDs)
                    leLIDs = msh.faceElems(1, faceIDs);
                    reLIDs = msh.faceElems(2, faceIDs);
                    
                    [rll, cll, vll] = find(FM{lfn, rfn, f2f}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn, f2f}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn, f2f}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn, f2f}{2, 2});
                    vll = 0.5 * vll;
                    vlr = 0.5 * vlr;
                    vrl = 0.5 * vrl; 
                    vrr = 0.5 * vrr;                    
                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v = msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr];
                    vx(s : e) = reshape(msh.faceNormalx(:, faceIDs) .* v, [l, 1]); 
                    vy(s : e) = reshape(msh.faceNormaly(:, faceIDs) .* v, [l, 1]);  
                    vz(s : e) = reshape(msh.faceNormalz(:, faceIDs) .* v, [l, 1]);
                end
                i = i + 1;
            end  
            Ax = sparse(r(1 : e), c(1 : e), vx(1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(1 : e), c(1 : e), vy(1 : e), msh.nLElems * nb, msh.nLElems * nb); 
            Az = sparse(r(1 : e), c(1 : e), vz(1 : e), msh.nLElems * nb, msh.nLElems * nb);            
        otherwise
            error('Wrong flux')
    end
else
    switch flux
        case 'l'
            faceIDs = [md.intLFaces{:, 3, :}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r; vz = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                f2f       = msh.faceR2L(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                normalz   = msh.faceNormalz(:, ID);
                J         = msh.faceJac(:, ID);
                
                leLID = msh.elemLID(faceElems(1)); lfn = faceNums(1);
                reLID = msh.elemLID(faceElems(2)); rfn = faceNums(2);
                ip    = normalx * wind(1) + normaly * wind(2) + normalz * wind(3);                 
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (leLID - 1) * nb];
                    v = J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{2, 1}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (reLID - 1) * nb, ce + (reLID - 1) * nb]; 
                    v = J * [FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 2}(:)];                  
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (leLID - 1) * nb, re + (reLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (reLID - 1) * nb, ce + (leLID - 1) * nb, ce + (reLID - 1) * nb];
                    v = 0.5 * J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 1}(:), FM{lfn, rfn, f2f}{2, 2}(:)];                   
                end
                vx(:, s : e) = normalx * v;
                vy(:, s : e) = normaly * v;
                vz(:, s : e) = normalz * v;                
            end
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);  
            Az = sparse(r(:, 1 : e), c(:, 1 : e), vz(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
        case 'r'
            faceIDs = [md.intLFaces{:, 3, :}];
            r = zeros(nb2, 3 * md.nIntLFaces); c = r; vx = r; vy = r; vz = r;
            for i = 1 : md.nIntLFaces
                ID        = faceIDs(i);
                faceElems = msh.faceElems(:, ID);
                faceNums  = msh.faceNums(:, ID);
                f2f       = msh.faceR2L(:, ID);
                normalx   = msh.faceNormalx(:, ID);
                normaly   = msh.faceNormaly(:, ID);
                normalz   = msh.faceNormalz(:, ID);
                J         = msh.faceJac(:, ID);
                
                leLID = msh.elemLID(faceElems(1)); lfn = faceNums(1);
                reLID = msh.elemLID(faceElems(2)); rfn = faceNums(2);
                ip = normalx * wind(1) + normaly * wind(2) + normalz * wind(3);
                if (ip > 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (reLID - 1) * nb, ce + (reLID - 1) * nb]; 
                    v = J * [FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 2}(:)];
                elseif (ip < 0)
                    s = e + 1; e = e + 2;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (leLID - 1) * nb];
                    v = J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{2, 1}(:)];
                else
                    s = e + 1; e = e + 4;
                    r(:, s : e) = [re + (leLID - 1) * nb, re + (leLID - 1) * nb, re + (reLID - 1) * nb, re + (reLID - 1) * nb];
                    c(:, s : e) = [ce + (leLID - 1) * nb, ce + (reLID - 1) * nb, ce + (leLID - 1) * nb, ce + (reLID - 1) * nb];
                    v = 0.5 * J * [FM{lfn, rfn, f2f}{1, 1}(:), FM{lfn, rfn, f2f}{1, 2}(:), FM{lfn, rfn, f2f}{2, 1}(:), FM{lfn, rfn, f2f}{2, 2}(:)];
                end
                vx(:, s : e) = normalx * v;
                vy(:, s : e) = normaly * v;
                vz(:, s : e) = normalz * v;                 
            end  
            Ax = sparse(r(:, 1 : e), c(:, 1 : e), vx(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(:, 1 : e), c(:, 1 : e), vy(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Az = sparse(r(:, 1 : e), c(:, 1 : e), vz(:, 1 : e), msh.nLElems * nb, msh.nLElems * nb);
        case 'c'  
            r = zeros(4 * md.nIntLFaces * nb2, 1); c = r; vx = r; vy = r; vz = r;   
            i = 1;
            while (i <= 2592)
                f2f = ceil(i / 324);
                if (md.nifs(f2f) == 0)
                    i = i + 324;
                    continue;
                end
                j       = i - 324 * (f2f - 1);
                lfn     = md.intLFaces{j, 1, f2f};
                rfn     = md.intLFaces{j, 2, f2f};
                faceIDs = md.intLFaces{j, 3, f2f};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    
                    [rll, cll, vll] = find(FM{lfn, rfn, f2f}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn, f2f}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn, f2f}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn, f2f}{2, 2});
                    vll = 0.5 * vll;
                    vlr = 0.5 * vlr;
                    vrl = 0.5 * vrl; 
                    vrr = 0.5 * vrr;                       
                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v = msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr];
                    vx(s : e) = reshape(msh.faceNormalx(:, faceIDs) .* v, [l, 1]); 
                    vy(s : e) = reshape(msh.faceNormaly(:, faceIDs) .* v, [l, 1]);  
                    vz(s : e) = reshape(msh.faceNormalz(:, faceIDs) .* v, [l, 1]);
                end
                i = i + 1;
            end    
            Ax = sparse(r(1 : e), c(1 : e), vx(1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Ay = sparse(r(1 : e), c(1 : e), vy(1 : e), msh.nLElems * nb, msh.nLElems * nb);
            Az = sparse(r(1 : e), c(1 : e), vz(1 : e), msh.nLElems * nb, msh.nLElems * nb); 
        otherwise
            error('Wrong flux')
    end
end


end








