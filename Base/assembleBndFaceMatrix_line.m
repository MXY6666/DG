% msh      : mesh of line element in 1D
% flux     : numerical flux
% bs       : basis function set data (a struct) or just the degree of 
%            polynomial (an integer)
% FM       : face matrix at reference geometry
% bndTypes : boundary types
% mtol     : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function A = assembleBndFaceMatrix_line(msh, flux, bs, md, FM, bndTypes, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 101) && (msh.type ~= 102)
    error('Wrong mesh type')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_line(msh);
end 

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference line
if (nargin < 5) || isempty(FM)
    FM = cell(2, 2);
    FM{1, 1} = computeFaceMatrix_refLine(1, 0, 'l', bs, mtol);
    FM{2, 2} = computeFaceMatrix_refLine(2, 0, 'l', bs, mtol);
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
faceIDs = [md.bndLFaces{1, bndLoc}, md.bndLFaces{2, bndLoc}];
nbf = length(faceIDs);
if (nbf == 0)
    A = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
else
    r = zeros(nbf * nb2, 1); c = r; v = r;
    e = 0;
    for i = 1 : nbf
        leLID = msh.faceElems(1, faceIDs(i));
        lfn   = msh.faceNums( 1, faceIDs(i));
        if (msh.type == 102)
            leLID = msh.elemLID(leLID);
        end
        
        [rll, cll, vll] = find(FM{lfn, lfn}{1, 1});
        if (flux == 'c')
            vll = 0.5 * vll;
        end        

        s = e + 1; e = e + length(vll);
        r(s : e) = rll + (leLID - 1) * nb;
        c(s : e) = cll + (leLID - 1) * nb;
        if (lfn == 1)
            v(s : e) = -vll;
        elseif (lfn == 2)
            v(s : e) = vll;
        else
            error('Wrong face number')
        end            
    end 
    A = sparse(r(1 : e), c(1 : e), v(1 : e), msh.nLElems * nb, msh.nLElems * nb);
end

end


