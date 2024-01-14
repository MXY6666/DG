% msh    : mesh 
% U      : coefficients of numerical solutions
% bs     : basis function set data
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function TVDmax = computeTVD(msh, U, bs, layout)

if (nargin < 3)
    error('Not enough arguments')
end

if (size(U, 1) ~= bs.nb || mod(size(U, 2), msh.nLElems) ~= 0)
    error('Wrong argument U')
end
  
if (nargin < 4) || isempty(layout)
    layout = 1;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

% Compute solution at all quadrature points 
uh = bs.phi{1} * U;

% Compute TVD in each element
nv = size(U, 2) / msh.nLElems;
TVD = sum(abs(uh(2 : end, :) - uh(1 : end - 1, :)));
if (layout == 1)
    TVD = TVD((0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
elseif (layout == 2)
    TVD = reshape(TVD, [nv, msh.nLElems]);
end

% Get all internal face IDs
isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.intLFaces, msh.bndLFaces{isPeriodic}];
else
    IndFaceIDs = msh.intLFaces;
end

% Compute maximum value of TVD in each element and its connected elements
nf    = length(IndFaceIDs);
leIDs = msh.faceElems(1, IndFaceIDs);
reIDs = msh.faceElems(2, IndFaceIDs);

TVDmax = TVD;
for i = 1 : nf
    leID = leIDs(i);
    reID = reIDs(i);
    
    TVDmax(:, leID) = max(TVDmax(:, leID), TVD(:, reID));
    TVDmax(:, reID) = max(TVDmax(:, reID), TVD(:, leID));
end

end


