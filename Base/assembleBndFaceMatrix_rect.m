% msh      : mesh of rectangular element in 2D
% flux     : numerical flux
% bs       : basis function set data (a struct) or just the degree of 
%            polynomial (an integer)
% FM       : face matrix at reference geometry
% bndTypes : boundary types
% mtol     : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function [Ax, Ay] = assembleBndFaceMatrix_rect(msh, flux, bs, md, FM, bndTypes, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 201) && (msh.type ~= 202)
    error('Wrong mesh type')
end

if (flux ~= 'l') && (flux ~= 'c')
    error('Wrong flux for the boundary face')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_rect(msh);
end 

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference square
if (nargin < 5) || isempty(FM)
    FM = cell(4, 4);
    for fn = 1 : 4
        FM{fn, fn} = computeFaceMatrix_refSquare(fn, fn, 1, 'l', bs, mtol);
    end
end

if (nargin < 6) || isempty(bndTypes)
    bndTypes = msh.bndTypes;
end

bndTypes = unique(bndTypes, 'stable');
[bndFound, bndLoc] = ismember(bndTypes, msh.bndTypes);
bndLoc = bndLoc(bndFound);

nb  = size(FM{1, 1}{1, 1}, 1); % number of basis functions
nb2 = nb * nb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble Ax 
nbf = 0;
for i = 1 : length(bndLoc)
    nbf = nbf + length(md.bndLFaces{1, bndLoc(i)}) + length(md.bndLFaces{2, bndLoc(i)});
end
if (nbf == 0)
    Ax = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
else
    r = zeros(nbf * nb2, 1); c = r; v = r;
    e = 0;
    for lfn = 1 : 2
        faceIDs = [md.bndLFaces{lfn, bndLoc}];
        if ~isempty(faceIDs)
            leLIDs = msh.faceElems(1, faceIDs);
            if (msh.type == 202)
                leLIDs = msh.elemLID(leLIDs);
            end
            J = msh.faceJac(:, faceIDs);
            
            [rll, cll, vll] = find(FM{lfn, lfn}{1, 1});
            if (flux == 'c')
                vll = 0.5 * vll;
            end            

            s = e + 1; l = length(faceIDs) * length(vll); e = e + l;
            r(s : e) = reshape(rll + (leLIDs - 1) * nb, [l, 1]);
            c(s : e) = reshape(cll + (leLIDs - 1) * nb, [l, 1]);
            if (lfn == 1)
                v(s : e) = reshape(-J .* vll, [l, 1]);
            else
                v(s : e) = reshape(J .* vll, [l, 1]);
            end          
        end          
    end  
    Ax = sparse(r(1 : e), c(1 : e), v(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble Ay
nbf = 0;
for i = 1 : length(bndLoc)
    nbf = nbf + length(md.bndLFaces{3, bndLoc(i)}) + length(md.bndLFaces{4, bndLoc(i)});
end
if (nbf == 0)
    Ay = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
else
    r = zeros(nbf * nb2, 1); c = r; v = r;
    e = 0;
    for lfn = 3 : 4
        faceIDs = [md.bndLFaces{lfn, bndLoc}];
        if ~isempty(faceIDs)
            leLIDs = msh.faceElems(1, faceIDs);
            if (msh.type == 202)
                leLIDs = msh.elemLID(leLIDs);
            end
            J = msh.faceJac(:, faceIDs);
            
            [rll, cll, vll] = find(FM{lfn, lfn}{1, 1});
            if (flux == 'c')
                vll = 0.5 * vll;
            end            

            s = e + 1; l = length(faceIDs) * length(vll); e = e + l;
            r(s : e) = reshape(rll + (leLIDs - 1) * nb, [l, 1]);
            c(s : e) = reshape(cll + (leLIDs - 1) * nb, [l, 1]);
            if (lfn == 3)
                v(s : e) = reshape(-J .* vll, [l, 1]);
            else
                v(s : e) = reshape(J .* vll, [l, 1]);
            end          
        end          
    end 
    Ay = sparse(r(1 : e), c(1 : e), v(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

end

