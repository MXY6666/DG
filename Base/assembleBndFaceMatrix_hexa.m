% msh      : mesh of linear hexahedral element in 3D
% flux     : numerical flux
% bs       : basis function set data (a struct) or just the degree of 
%            polynomial (an integer)
% FM       : face matrix at reference geometry
% bndTypes : boundary types
% mtol     : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function [Ax, Ay, Az] = assembleBndFaceMatrix_hexa(msh, flux, bs, md, FM, bndTypes, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 311) && (msh.type ~= 312)
    error('Wrong mesh type')
end

if (flux ~= 'l') && (flux ~= 'c')
    error('Wrong flux for the boundary face')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_hexa(msh);
end 

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference cube
if (nargin < 5) || isempty(FM)
    FM = cell(6, 6, 1);
    for fn = 1 : 6
        FM{fn, fn, 1} = computeFaceMatrix_refCube(fn, fn, 1, 'l', bs, mtol);
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
% Assemble
nbf = sum(msh.nBndLFaces(bndLoc));
if (nbf == 0)
    Ax = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
    Ay = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
    Az = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
else
    r = zeros(nbf * nb2, 1); c = r; vx = r; vy = r; vz = r;
    e = 0;
    for lfn = 1 : 6
        faceIDs = [md.bndLFaces{lfn, bndLoc}];
        if ~isempty(faceIDs)
            leLIDs = msh.faceElems(1, faceIDs);
            if (msh.type == 312)
                leLIDs = msh.elemLID(leLIDs);
            end
            J = msh.faceJac(:, faceIDs);
            
            [rll, cll, vll] = find(FM{lfn, lfn}{1, 1});
            if (flux == 'c')
                vll = 0.5 * vll;
            end            
            
            s = e + 1; l = length(faceIDs) * length(vll); e = e + l;
            r(s : e)  = reshape(rll + (leLIDs - 1) * nb, [l, 1]);
            c(s : e)  = reshape(cll + (leLIDs - 1) * nb, [l, 1]);
            vx(s : e) = reshape(J .* msh.faceNormalx(:, faceIDs) .* vll, [l, 1]);
            vy(s : e) = reshape(J .* msh.faceNormaly(:, faceIDs) .* vll, [l, 1]);
            vz(s : e) = reshape(J .* msh.faceNormalz(:, faceIDs) .* vll, [l, 1]);
        end
    end
    Ax = sparse(r(1 : e), c(1 : e), vx(1 : e), msh.nLElems * nb, msh.nLElems * nb);
    Ay = sparse(r(1 : e), c(1 : e), vy(1 : e), msh.nLElems * nb, msh.nLElems * nb);
    Az = sparse(r(1 : e), c(1 : e), vz(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

end


