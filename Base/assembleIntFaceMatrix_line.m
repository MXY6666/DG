% msh  : mesh of line element in 1D
% flux : numerical flux
% bs   : basis function set data (a struct) or just the degree of 
%        polynomial (an integer)
% md   : mesh data
% FM   : face matrix at reference geometry
% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function A = assembleIntFaceMatrix_line(msh, flux, bs, md, FM, mtol)

if (nargin < 3)
    error('Not enough arguments')
end   

if (msh.type ~= 101) && (msh.type ~= 102)
    error('Wrong mesh type')
end

if (nargin < 4) || isempty(md)
    md = computeMeshData_line(msh);
end 

if (nargin < 6) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute face matrix at reference line
if (nargin < 5) || isempty(FM)
    FM = cell(2, 2);
    FM{2, 1} = computeFaceMatrix_refLine(2, 1, flux, bs, mtol);
end

switch flux
    case {'l', 'c'}
        nb = size(FM{2, 1}{1, 1}, 1);
    case 'r'
        nb = size(FM{2, 1}{1, 2}, 1);
    otherwise
        error('Wrong flux') 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble
faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
if (msh.type == 102)
    leLIDs = msh.elemLID(leLIDs);
    reLIDs = msh.elemLID(reLIDs);
end
switch (flux)
    case 'l'
        [rll, cll, vll] = find(FM{2, 1}{1, 1});
        [rrl, crl, vrl] = find(FM{2, 1}{2, 1});
        
        r = [rll + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb];
        c = [cll + (leLIDs - 1) * nb; crl + (leLIDs - 1) * nb];
        v = repmat([vll; vrl], [1, length(faceIDs)]);        
    case 'r'
        [rlr, clr, vlr] = find(FM{2, 1}{1, 2});
        [rrr, crr, vrr] = find(FM{2, 1}{2, 2});
        
        r = [rlr + (leLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
        c = [clr + (reLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
        v = repmat([vlr; vrr], [1, length(faceIDs)]);         
    case 'c'
        [rll, cll, vll] = find(FM{2, 1}{1, 1});
        [rlr, clr, vlr] = find(FM{2, 1}{1, 2});
        [rrl, crl, vrl] = find(FM{2, 1}{2, 1});
        [rrr, crr, vrr] = find(FM{2, 1}{2, 2});
        
        r = [rll + (leLIDs - 1) * nb; rlr + (leLIDs - 1) * nb; rrl + (reLIDs - 1) * nb; rrr + (reLIDs - 1) * nb];
        c = [cll + (leLIDs - 1) * nb; clr + (reLIDs - 1) * nb; crl + (leLIDs - 1) * nb; crr + (reLIDs - 1) * nb];
        v = repmat([vll; vlr; vrl; vrr], [1, length(faceIDs)]);               
    otherwise
        error('Wrong flux')
end

A = sparse(r, c, v, msh.nLElems * nb, msh.nLElems * nb);

end



