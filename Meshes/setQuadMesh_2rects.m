% The computational domain is a rectangle (described with rect1) with a 
% small rectangle (described with rect2) embedded in it.
% Copyright (C) Xiangyi Meng.

function msh = setQuadMesh_2rects(rect1, rect2, N, bcs, type, maxLevel, quad, bs)

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
    type = 211;
end
if (type ~= 211) && (type ~= 212)
    error('wrong mesh type')
end

if (nargin < 6)  || isempty(maxLevel) || (type == 211) 
    maxLevel = 0;
end
if (type == 212) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 7) || isempty(quad)
    quad = GaussQuadratureRule_quad([2, 2], 201);
end

if (nargin < 8) || isempty(bs)
    quad1 = GaussQuadratureRule_line(2, 101);
    bs = setBasisFunctionSet_square(quad1, quad, 201, 2);
end

msh = setQuadMesh_rect(rect1, N, bcs(1 : 4), type, maxLevel, [0, 0], quad, bs);
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

% Store the IDs of nodes which will be deleted
dnIDs = (Nx1 + 2 : Nx1 + Nx2)' + (Ny1 + 1 : Ny1 + Ny2 - 1) * (N(1) + 1);
dnIDs = dnIDs(:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change the information of the internal solid wall faces
% the left internal solid wall faces
IDs = msh.elemFaces(1, deIDs(1, :));
msh.faceType(IDs)     = bcs(5);
msh.faceElems(2, IDs) = 0;
msh.faceNums(2, IDs)  = 0;
% the right internal solid wall faces
IDs = msh.elemFaces(2, deIDs(Nx2, :));
msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
msh.faceType(IDs)     = bcs(5);
temp = msh.faceNodes(1, IDs);
msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
msh.faceNodes(2, IDs) = temp;  
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
msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
msh.faceType(IDs)     = bcs(5);
temp = msh.faceNodes(1, IDs);
msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
msh.faceNodes(2, IDs) = temp;  
msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
msh.faceElems(2, IDs) = 0;
msh.faceNums(1, IDs)  = 3;
msh.faceNums(2, IDs)  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refresh element-node, element-face and face-element and face-node connections
deIDs = deIDs(:)';
% Re-index the remaining elements, faces and nodes
reIDs = 1 : msh.nElems;
rfIDs = 1 : msh.nFaces;
rnIDs = 1 : msh.nNodes;
reIDs(deIDs) = []; % remaining element old IDs
rfIDs(dfIDs) = []; % remaining face old IDs
rnIDs(dnIDs) = []; % remaining node old IDs
nre = length(reIDs);
nrf = length(rfIDs);
nrn = length(rnIDs);

neIDs = zeros(1, msh.nElems); 
nfIDs = zeros(1, msh.nFaces); 
nnIDs = zeros(1, msh.nNodes);
neIDs(reIDs) = 1 : nre; % new element IDs
nfIDs(rfIDs) = 1 : nrf; % new face IDs
nnIDs(rnIDs) = 1 : nrn; % new node IDs

% Refresh element-node and element-face connections
msh.elemNodes(:, reIDs) = nnIDs(msh.elemNodes(:, reIDs));
msh.elemFaces(:, reIDs) = nfIDs(msh.elemFaces(:, reIDs));

% Refresh face-element and face-node connections
msh.faceNodes(:, rfIDs) = nnIDs(msh.faceNodes(:, rfIDs));
msh.faceElems(1, rfIDs) = neIDs(msh.faceElems(1, rfIDs));
rfIDs_int = rfIDs(msh.faceType(rfIDs) == 0 | msh.faceType(rfIDs) == 1);
msh.faceElems(2, rfIDs_int) = neIDs(msh.faceElems(2, rfIDs_int));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally remove the deleted elements, faces and nodes from the mesh
msh.elemSize(deIDs)       = [];
msh.elemDiam(deIDs)       = [];
msh.elemNodes(:, deIDs)   = [];
msh.elemFaces(:, deIDs)   = [];
msh.elemGPx(:, deIDs)     = []; 
msh.elemGPy(:, deIDs)     = [];
msh.elemJac(:, deIDs)     = [];
msh.elemJxix(:, deIDs)    = []; 
msh.elemJxiy(:, deIDs)    = [];
msh.elemJetax(:, deIDs)   = []; 
msh.elemJetay(:, deIDs)   = [];
msh.elemIsCurved(deIDs)   = [];
msh.faceNormalx(:, dfIDs) = [];
msh.faceNormaly(:, dfIDs) = [];
msh.faceSize(dfIDs)       = [];
msh.faceType(dfIDs)       = [];
msh.faceNodes(:, dfIDs)   = [];
msh.faceElems(:, dfIDs)   = [];
msh.faceNums(:, dfIDs)    = [];
msh.faceJac(:, dfIDs)     = [];
msh.faceIsCurved(dfIDs)   = [];
msh.nodeCoor(:, dnIDs)    = [];
msh.massMatInv(deIDs)     = [];
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
% Update the IDs of leaf elements and faces
if (type == 212)
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




