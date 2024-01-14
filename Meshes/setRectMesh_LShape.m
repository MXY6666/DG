% domain     : a 4D vector 
% N          : a 2D vector
% des        : a 3D vector, the fist and second entries store the length of the 
%              deleted rectangle, the third entry of 1(left-bottom), 2(right-bottom), 
%              3(left-top) or 4(right-top) represents the location of the deleted rectangle.
% bcs        : 1 * 4 or 4 * 1 vector containing the four boundary types
% type = 201 : ignore the adaptive information      
% type = 202 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% Copyright (C) Xiangyi Meng.

function msh = setRectMesh_LShape(domain, N, des, bcs, type, maxLevel)

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

msh = setRectMesh_rect(domain, N, bcs, type, maxLevel);
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
        
        % Change the information of some faces
        IDs = msh.elemFaces(2, deIDs(end, :));
        msh.faceNormalx(IDs)  = -1;
        msh.faceType(IDs)     = bcs(3);        
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 1;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(4, deIDs(:, end));
        msh.faceNormaly(IDs)  = -1;
        msh.faceType(IDs)     = bcs(3);         
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
        
        % Change the information of some faces
        IDs = msh.elemFaces(1, deIDs(1, :));
        msh.faceType(IDs)     = bcs(3);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;
        
        IDs = msh.elemFaces(4, deIDs(:, end));
        msh.faceNormaly(IDs)  = -1;
        msh.faceType(IDs)     = bcs(3);        
        msh.faceElems(1, IDs) = msh.faceElems(2, IDs);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 1, IDs) = 3;
        msh.faceNums( 2, IDs) = 0;      
    case 3
        % Store the IDs of elements which will be deleted
        deIDs = ((N(2) - Ny) * N(1) + 1 : (N(2) - Ny) * N(1) + Nx)' + (0 : Ny - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(1, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(4, deIDs);   
        
        % Change the information of some faces
        IDs = msh.elemFaces(2, deIDs(end, :));
        msh.faceNormalx(IDs)  = -1;
        msh.faceType(IDs)     = bcs(4);        
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
        deIDs = ((N(2) - Ny + 1) * N(1) - Nx + 1 : (N(2) - Ny + 1) * N(1))' + (0 : Ny - 1) * N(1);
        
        % Store the IDs of faces which will be deleted
        dfIDs(1 : nde)       = msh.elemFaces(2, deIDs);
        dfIDs(nde + 1 : end) = msh.elemFaces(4, deIDs);    
        
        % Change the information of some faces
        IDs = msh.elemFaces(1, deIDs(1, :));
        msh.faceType(IDs) = bcs(4);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;

        IDs = msh.elemFaces(3, deIDs(:, 1));
        msh.faceType(IDs) = bcs(4);
        msh.faceElems(2, IDs) = 0;
        msh.faceNums( 2, IDs) = 0;      
    otherwise
        error('Wrong argument des')
end

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


