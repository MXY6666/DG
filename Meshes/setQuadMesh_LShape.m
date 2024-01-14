% domain     : a 4D vector 
% N          : a 2D vector
% des        : a 3D vector, the fist and second entries store the length of the 
%              deleted rectangle, the third entry of 1(left-bottom), 2(right-bottom), 
%              3(left-top) or 4(right-top) represents the location of the deleted rectangle.
% bcs        : 1 * 4 or 4 * 1 vector containing the four boundary types
% type = 211 : ignore the adaptive information        
% type = 212 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% quad       : Gauss quadrature rule for a quadrilateral
% bs         : basis function set
% Copyright (C) Xiangyi Meng.

function msh = setQuadMesh_LShape(domain, N, des, bcs, type, maxLevel, quad, bs)

if (nargin < 3)
    error('Not enough arguments')
end

if (length(domain) ~= 4)
    error('Wrong size of argument domain')
end

if (length(N) ~= 2)
    error('Wrong size of argument N')
end

if (length(des) ~= 3)
    error('Wrong size of argument des')
end

if (domain(2) <= domain(1)) || (domain(4) <= domain(3))
    error('Wrong argument domain')
end

if (des(3) ~= 1 && des(3) ~= 2 && des(3) ~= 3 && des(3) ~= 4)
    error('The third entry of argument des can only be 1, 2, 3 or 4')
end

fracx = des(1) / (domain(2) - domain(1));
fracy = des(2) / (domain(4) - domain(3));    
Nx = round(N(1) * fracx);
Ny = round(N(2) * fracy);
if (abs(N(1) * fracx - Nx) > 1e-12 || abs(N(2) * fracy - Ny) > 1e-12)
    error('This number of elements will cause hanging node in the corner')
end

if (nargin < 4)  || isempty(bcs)
    bcs = [6, 7, 8, 8];
end
if (length(bcs) ~= 4)
    error('Wrong boundary conditions')
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

msh = setQuadMesh_rect(domain, N, bcs, type, maxLevel, [0, 0], quad, bs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nde   = Nx * Ny; % number of deleted elements
ndf   = 2 * nde; % number of deleted faces
dfIDs = zeros(1, ndf); % face IDs which will be deleted
switch des(3)
    case 1
        % Store the IDs of elements which will be deleted
        deIDs = (1 : Nx)' + (0 : Ny - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(1, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(3, deIDs);
        
        % Store the IDs of nodes which will be deleted
        dnIDs = (1 : Nx)' + (0 : Ny - 1) * (N(1) + 1);
        
        % Change the information of some faces
        IDs = msh.elemFaces(2, deIDs(end, :));
        msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
        msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
        msh.faceType(IDs)     = bcs(3);
        temp = msh.faceNodes(1, IDs);
        msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
        msh.faceNodes(2, IDs) = temp;        
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 1;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(4, deIDs(:, end));
        msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
        msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
        msh.faceType(IDs)     = bcs(3);  
        temp = msh.faceNodes(1, IDs);
        msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
        msh.faceNodes(2, IDs) = temp;         
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 3;
        msh.faceNums( 2, IDs) = 0;   
    case 2
        % Store the IDs of elements which will be deleted
        deIDs = (N(1) - Nx + 1 : N(1))' + (0 : Ny - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(2, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(3, deIDs);   
        
        % Store the IDs of nodes which will be deleted
        dnIDs = (N(1) - Nx + 2 : N(1) + 1)' + (0 : Ny - 1) * (N(1) + 1);        
        
        % Change the information of some faces
        IDs = msh.elemFaces(1, deIDs(1, :));
        msh.faceType(IDs)     = bcs(3);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(4, deIDs(:, end));
        msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
        msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
        msh.faceType(IDs)     = bcs(3);  
        temp = msh.faceNodes(1, IDs);
        msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
        msh.faceNodes(2, IDs) = temp;         
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 3;
        msh.faceNums( 2, IDs) = 0;        
    case 3
        % Store the IDs of elements which will be deleted
        deIDs = (1 : Nx)' + (N(2) - Ny : N(2) - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(1, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(4, deIDs);   
        
        % Store the IDs of nodes which will be deleted
        dnIDs = (1 : Nx)' + (N(2) - Ny + 1 : N(2)) * (N(1) + 1);        
        
        % Change the information of some faces
        IDs = msh.elemFaces(2, deIDs(end, :));
        msh.faceNormalx(IDs)  = -msh.faceNormalx(IDs);
        msh.faceNormaly(IDs)  = -msh.faceNormaly(IDs);
        msh.faceType(IDs)     = bcs(4);
        temp = msh.faceNodes(1, IDs);
        msh.faceNodes(1, IDs) = msh.faceNodes(2, IDs);
        msh.faceNodes(2, IDs) = temp;        
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 1;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(3, deIDs(:, 1));
        msh.faceType(IDs)     = bcs(4);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;
    case 4
        % Store the IDs of elements which will be deleted
        deIDs = (N(1) - Nx + 1 : N(1))' + (N(2) - Ny : N(2) - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(2, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(4, deIDs);
        
        % Store the IDs of nodes which will be deleted
        dnIDs = (N(1) - Nx + 2 : N(1) + 1)' + (N(2) - Ny + 1 : N(2)) * (N(1) + 1);         
        
        % Change the information of some faces
        IDs = msh.elemFaces(1, deIDs(1, :));
        msh.faceType(IDs)     = bcs(4);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(3, deIDs(:, 1));
        msh.faceType(IDs)     = bcs(4);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;     
    otherwise
        error('Wrong argument des')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refresh element-node, element-face and face-element and face-node connections
deIDs = deIDs(:)';
dnIDs = dnIDs(:)';
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
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = find(msh.faceType == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Augment information for adaptive mesh
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


