% msh  : mesh of rectangular element in 2D
% flux : numerical flux
% bs   : basis function set data (a struct) or just the degree of 
%        polynomial (an integer)
% md   : mesh data
% FM   : face matrix at reference geometry
% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function [Ax, Ay] = assembleIntFaceMatrix_rect(msh, flux, bs, md, FM, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 201) && (msh.type ~= 202)
    error('Wrong mesh type')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_rect(msh);
end 

if (nargin < 6) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference square
if (nargin < 5) || isempty(FM)
    if (msh.type == 201)
        FM = cell(4, 4);
    else
        FM = cell(12, 12);
    end
    for fk = 1 : 2
        [lfn, rfn] = getFaceNumber_rect(fk);
        FM{lfn, rfn} = computeFaceMatrix_refSquare(lfn, rfn, 1, flux, bs, mtol);
    end    
    if (msh.type == 202)
        for fk = 3 : 10
            [lfn, rfn] = getFaceNumber_rect(fk);
            FM{lfn, rfn} = computeFaceMatrix_refSquare(lfn, rfn, 1, flux, bs, mtol);
        end
    end
end

switch flux
    case {'l', 'c'}
        nb = size(FM{2, 1}{1, 1}, 1);
    case 'r'
        nb = size(FM{2, 1}{1, 2}, 1);
    otherwise
        error('Wrong flux') 
end
nb2 = nb * nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble Ax first 
if (msh.type == 201)
    faceIDs = md.intLFaces{1, 3};
    leLIDs  = msh.faceElems(1, faceIDs);
    reLIDs  = msh.faceElems(2, faceIDs);
    switch flux
        case 'l'
            [rll, cll, vll] = find(FM{2, 1}{1, 1});
            [rrl, crl, vrl] = find(FM{2, 1}{2, 1});  
            
            r = [rll + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb];
            c = [cll + (leLIDs - 1) * nb; crl + (leLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vll; vrl];        
        case 'r'
            [rlr, clr, vlr] = find(FM{2, 1}{1, 2});
            [rrr, crr, vrr] = find(FM{2, 1}{2, 2});
            
            r = [rlr + (leLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
            c = [clr + (reLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vlr; vrr];            
        case 'c'
            [rll, cll, vll] = find(FM{2, 1}{1, 1});
            [rlr, clr, vlr] = find(FM{2, 1}{1, 2});
            [rrl, crl, vrl] = find(FM{2, 1}{2, 1});
            [rrr, crr, vrr] = find(FM{2, 1}{2, 2});    
            
            r = [rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
            c = [cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr];
        otherwise
            error('Wrong flux')
    end 
    Ax = sparse(r, c, v, msh.nLElems * nb, msh.nLElems * nb);
else
    e = 0;
    switch flux
        case 'l'
            r = zeros(2 * md.nIntLFacesx * nb2, 1); c = r; v = r;           
            for i = [1, 3 : 6]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs)); 
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1}); 

                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vrl)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; crl + (leLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vll; vrl], [l, 1]);
                end
            end
        case 'r'
            r = zeros(2 * md.nIntLFacesx * nb2, 1); c = r; v = r;
            for i = [1, 3 : 6]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2});
                
                    s = e + 1; l = length(faceIDs) * (length(vlr) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rlr + (leLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([clr + (reLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vlr; vrr], [l, 1]);
                end
            end            
        case 'c'
            r = zeros(4 * md.nIntLFacesx * nb2, 1); c = r; v = r;      
            for i = [1, 3 : 6]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2}); 
                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr], [l, 1]);
                end
            end
        otherwise
            error('Wong flux')
    end
    Ax = sparse(r(1 : e), c(1 : e), v(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble Ay 
if (msh.type == 201)
    faceIDs = md.intLFaces{2, 3};
    leLIDs  = msh.faceElems(1, faceIDs);
    reLIDs  = msh.faceElems(2, faceIDs);
    switch flux
        case 'l'
            [rll, cll, vll] = find(FM{4, 3}{1, 1});
            [rrl, crl, vrl] = find(FM{4, 3}{2, 1});
            
            r = [rll + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb];
            c = [cll + (leLIDs - 1) * nb; crl + (leLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vll; vrl];         
        case 'r'
            [rlr, clr, vlr] = find(FM{4, 3}{1, 2});
            [rrr, crr, vrr] = find(FM{4, 3}{2, 2});
            
            r = [rlr + (leLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
            c = [clr + (reLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vlr; vrr];            
        case 'c'
            [rll, cll, vll] = find(FM{4, 3}{1, 1});
            [rlr, clr, vlr] = find(FM{4, 3}{1, 2});
            [rrl, crl, vrl] = find(FM{4, 3}{2, 1});
            [rrr, crr, vrr] = find(FM{4, 3}{2, 2});
            
            r = [rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
            c = [cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
            v = msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr];              
        otherwise
            error('Wrong flux')
    end
    Ay = sparse(r, c, v, msh.nLElems * nb, msh.nLElems * nb);
else
    e = 0;
    switch flux
        case 'l'
            r = zeros(2 * md.nIntLFacesy * nb2, 1); c = r; v = r;           
            for i = [2, 7 : 10]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs)); 
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1}); 

                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vrl)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; crl + (leLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vll; vrl], [l, 1]);
                end
            end
        case 'r'
            r = zeros(2 * md.nIntLFacesy * nb2, 1); c = r; v = r;
            for i = [2, 7 : 10]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2});
                
                    s = e + 1; l = length(faceIDs) * (length(vlr) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rlr + (leLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([clr + (reLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vlr; vrr], [l, 1]);
                end
            end            
        case 'c'
            r = zeros(4 * md.nIntLFacesy * nb2, 1); c = r; v = r;      
            for i = [2, 7 : 10]
                lfn     = md.intLFaces{i, 1};
                rfn     = md.intLFaces{i, 2};
                faceIDs = md.intLFaces{i, 3};
                if ~isempty(faceIDs)
                    leLIDs = msh.elemLID(msh.faceElems(1, faceIDs));
                    reLIDs = msh.elemLID(msh.faceElems(2, faceIDs));
                    
                    [rll, cll, vll] = find(FM{lfn, rfn}{1, 1});
                    [rlr, clr, vlr] = find(FM{lfn, rfn}{1, 2});
                    [rrl, crl, vrl] = find(FM{lfn, rfn}{2, 1});
                    [rrr, crr, vrr] = find(FM{lfn, rfn}{2, 2}); 
                    
                    s = e + 1; l = length(faceIDs) * (length(vll) + length(vlr) + length(vrl) + length(vrr)); e = e + l;
                    r(s : e) = reshape([rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb], [l, 1]);
                    c(s : e) = reshape([cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb], [l, 1]);
                    v(s : e) = reshape(msh.faceJac(:, faceIDs) .* [vll; vlr; vrl; vrr], [l, 1]);
                end
            end
        otherwise
            error('Wong flux')
    end
    Ay = sparse(r(1 : e), c(1 : e), v(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

end
