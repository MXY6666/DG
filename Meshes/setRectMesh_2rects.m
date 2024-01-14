% The computational domain is a rectangle (described with rect1) with a 
% small rectangle (described with rect2) embedded in it.
% Copyright (C) Xiangyi Meng.

function msh = setRectMesh_2rects(rect1, rect2, N, bcs, type, maxLevel)

if (nargin < 3)
    error('Not enough arguments')
end

if (length(rect1) ~= 4)
    error('Wrong size of argument rect1')
end

if (length(rect2) ~= 4)
    error('Wrong size of argument rect2')
end

if (length(N) ~= 2)
    error('Wrong size of argument N')
end

if ~(rect1(1) < rect2(1) && rect2(1) < rect2(2) && rect2(2) < rect1(2))
    error('Wrong domain boundary')
end

if ~(rect1(3) < rect2(3) && rect2(3) < rect2(4) && rect2(4) < rect1(4))
    error('Wrong domain boundary')
end

frac1 = (rect2(1) - rect1(1)) / (rect1(2) - rect1(1));
frac2 = (rect2(2) - rect2(1)) / (rect1(2) - rect1(1));
frac3 = (rect1(2) - rect2(2)) / (rect1(2) - rect1(1));
Nx1 = round(N(1) * frac1);
Nx2 = round(N(1) * frac2);
Nx3 = round(N(1) * frac3);
if (abs(Nx1 - N(1) * frac1) > 1e-12 || abs(Nx2 - N(1) * frac2) > 1e-12 || abs(Nx3 - N(1) * frac3) > 1e-12)
    error('The initial mesh data may cause hanging node')
end

frac1 = (rect2(3) - rect1(3)) / (rect1(4) - rect1(3));
frac2 = (rect2(4) - rect2(3)) / (rect1(4) - rect1(3));
frac3 = (rect1(4) - rect2(4)) / (rect1(4) - rect1(3));
Ny1 = round(N(2) * frac1);
Ny2 = round(N(2) * frac2);
Ny3 = round(N(2) * frac3);
if (abs(Ny1 - N(2) * frac1) > 1e-12 || abs(Ny2 - N(2) * frac2) > 1e-12 || abs(Ny3 - N(2) * frac3) > 1e-12)
    error('The initial mesh data may cause hanging node')
end

if (nargin < 4)  || isempty(bcs)
    bcs = [6, 7, 8, 8, 8];
end
if (length(bcs) ~= 5)
    error('Wrong boundary conditions')
end
if (bcs(5) == 0)
    error('The internal solid wall boundary cannot be indexed periodic')
end

if (nargin < 5)  || isempty(type)
    type = 201;
end
if (type ~= 201) && (type ~= 202) 
    error('wrong mesh type')
end

if (nargin < 6)  || isempty(maxLevel) || (type == 201) 
    maxLevel = 0;
end
if (type == 202) && (maxLevel == 0)
    maxLevel = 3;
end

msh = setRectMesh_rect(rect1, N, bcs(1 : 4), type, maxLevel);
msh.bndTypes = unique(bcs, 'stable');
msh.rect = rect2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store the IDs of elements which will be deleted
nde   = Nx2 * Ny2;
deIDs = (Nx1 + 1 : Nx1 + Nx2)' + (Ny1 : Ny1 + Ny2 - 1) * N(1);

% Store the IDs of faces which will be deleted
ndf   = 2 * nde - Nx2 - Ny2;
dfIDs = zeros(1, ndf);
dfIDs(1 : nde - Ny2)       = msh.elemFaces(2, deIDs(1 : Nx2 - 1, :)); % the vertical faces
dfIDs(nde - Ny2 + 1 : end) = msh.elemFaces(4, deIDs(:, 1 : Ny2 - 1)); % the horizontal faces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the information of the internal solid wall faces
% the left internal solid wall faces
IDs = msh.elemFaces(1, deIDs(1, :));
msh.faceType(IDs)     = bcs(5);
msh.faceElems(2, IDs) = 0;
msh.faceNums(2, IDs)  = 0;
% the right internal solid wall faces
IDs = msh.elemFaces(2, deIDs(Nx2, :));
msh.faceNormalx(IDs)  = -1;
msh.faceType(IDs)     = bcs(5);
msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
msh.faceElems(2, IDs) = 0;
msh.faceNums(1, IDs)  = 1;
msh.faceNums(2, IDs)  = 0;
% the bottom internal solid wall faces
IDs = msh.elemFaces(3, deIDs(:, 1));
msh.faceType(IDs)     = bcs(5);
msh.faceElems(2, IDs) = 0;
msh.faceNums(2, IDs)  = 0;
% the top internal solid wall faces
IDs = msh.elemFaces(4, deIDs(:, Ny2));
msh.faceNormaly(IDs)  = -1;
msh.faceType(IDs)     = bcs(5);
msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
msh.faceElems(2, IDs) = 0;
msh.faceNums(1, IDs)  = 3;
msh.faceNums(2, IDs)  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refresh element-face and face-element connections
deIDs = deIDs(:)';
% Re-index the remaining elements, faces 
reIDs = 1 : msh.nElems;
rfIDs = 1 : msh.nFaces;
reIDs(deIDs) = []; % remaining element old IDs
rfIDs(dfIDs) = []; % remaining face old IDs
nre = length(reIDs);
nrf = length(rfIDs);

neIDs = zeros(1, msh.nElems); 
nfIDs = zeros(1, msh.nFaces); 
neIDs(reIDs) = 1 : nre; % new element IDs
nfIDs(rfIDs) = 1 : nrf; % new face IDs

% Refresh element-face connections
msh.elemFaces(:, reIDs) = nfIDs(msh.elemFaces(:, reIDs));

% Refresh face-element connections
msh.faceElems(1, rfIDs) = neIDs(msh.faceElems(1, rfIDs));
rfIDs_int = rfIDs(msh.faceType(rfIDs) == 0 | msh.faceType(rfIDs) == 1);
msh.faceElems(2, rfIDs_int) = neIDs(msh.faceElems(2, rfIDs_int));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally remove the deleted elements, faces and nodes from the mesh
msh.elemCenter(:, deIDs)  = [];
msh.elemLength(:, deIDs)  = [];
msh.elemSize(deIDs)       = [];
msh.elemFaces(:, deIDs)   = [];
msh.elemJac(:, deIDs)     = [];
msh.elemJxix(:, deIDs)    = [];
msh.elemJetay(:, deIDs)   = [];
msh.faceNormalx(:, dfIDs) = [];
msh.faceNormaly(:, dfIDs) = [];
msh.faceSize(dfIDs)       = [];
msh.faceType(dfIDs)       = [];
msh.faceElems(:, dfIDs)   = [];
msh.faceNums(:, dfIDs)    = [];
msh.faceJac(:, dfIDs)     = [];
msh.nElems     = nre;
msh.nFaces     = nrf;
msh.nLElems    = nre;
msh.nLFaces    = nrf;
msh.LElems     = 1 : nre;
msh.LFaces     = 1 : nrf;
msh.intLFaces  = find(msh.faceType == 0);
msh.nIntLFaces = length(msh.intLFaces);
msh.bndLFaces  = cell(1, length(msh.bndTypes));
msh.nBndLFaces = zeros(1, length(msh.bndTypes));
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = find(msh.faceType == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Augment information for adaptive mesh
if (type == 202)
    msh.elemLevel    = zeros(1, nre);
    msh.elemLID      = 1 : nre;
    msh.elemParent   = zeros(1, nre);
    msh.elemChildren = zeros(4, nre);
    
    msh.faceLevel    = zeros(1, nrf);
    msh.faceLID      = 1 : nrf;
    msh.faceParent   = zeros(1, nrf);
    msh.faceChildren = zeros(2, nrf);
end

end




