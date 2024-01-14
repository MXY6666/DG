% List all the field of a mesh. Note, for a specific mesh type, some fileds 
% may be ignored
% face type 0  : internal
% face type 1  : periodic
% face type 2  : Dirichlet
% face type 3  : Neumann
% face type 4  : subsonic inflow
% face type 5  : subsonic outflow
% face type 6  : supersonic inflow
% face type 7  : supersonic outflow
% face type 8  : slip solid wall (symmetry or reflective wall)
% face type 9  : no-slip, adiabatic solid wall
% face type 10 : no-slip, isothermal solid wall
% face type 11 : exterior interpolation 
% Note, you may also mark any index for face type to satisfy your own
% needs.
% Copyright (C) Xiangyi Meng.

function msh = getEmptyMesh

% some general information of the mesh
msh.dm         = []; % domain information
msh.N          = []; % information of number of elements for initial struct mesh
msh.type       = []; % type information
msh.maxLevel   = 0;  % maximum tree level of the mesh
msh.bndTypes   = []; % all the types of boundaries
msh.nElems     = 0;  % number of elements
msh.nFaces     = 0;  % number of faces
msh.nEdges     = []; % number of edges
msh.nNodes     = []; % number of nodes
msh.nLElems    = 0;  % number of leaf elements
msh.nLFaces    = 0;  % number of leaf faces
msh.nIntLFaces = 0;  % number of internal leaf faces
msh.nBndLFaces = 0;  % number of boundary leaf faces
msh.LElems     = []; % leaf element IDs
msh.LFaces     = []; % leaf face IDs
msh.intLFaces  = []; % internal leaf face IDs
msh.bndLFaces  = []; % boundary leaf face IDs
msh.nGPs       = []; % number of Gauss points in each direction for each element
msh.massMatInv = []; % inverse of mass matrix in each element
%**************************************************************************
% metric information of elements
msh.elemCenter = []; % center 
msh.elemLength = []; % length in x, y, z direction for struct mesh
msh.elemSize   = []; % size(length, area or volume)
msh.elemDiam   = []; % diameter

% topology information of elements
msh.elemNodes = []; % node IDs of each element
msh.elemEdges = []; % edge IDs of each element
msh.elemFaces = []; % face IDs of each element

% Augment metric information of elements for adaptive mesh
msh.elemLevel  = []; % tree level
msh.elemLID    = []; % ID among the leaf elements

% Augment topology information of elements for adaptive mesh
msh.elemParent   = []; % parent ID of each element
msh.elemChildren = []; % child IDs of each element

% Store other useful information for computation
msh.elemGPx      = []; % x coordinate of Gauss points for each element
msh.elemGPy      = []; % y coordinate of Gauss points for each element
msh.elemGPz      = []; % z coordinate of Gauss points for each element
msh.elemJac      = []; % Jacobian of the transformation from the referce element to a physical element
msh.elemJxix     = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJxiy     = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJxiz     = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJetax    = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJetay    = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJetaz    = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJzetax   = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJzetay   = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemJzetaz   = []; % multiplication of Jacobian and the entry of inverse Jacobi matrix
msh.elemIsCurved = []; % An element is indexed curved if any of its face is indexed curved
%**************************************************************************
% metric information of faces
msh.faceNormalx = []; % x component of face normal vector
msh.faceNormaly = []; % y component of face normal vector
msh.faceNormalz = []; % z component of face normal vector
msh.faceSize    = []; % size
msh.faceType    = []; % type

% topology information of faces
msh.faceNodes = []; % node IDs of each face
msh.faceEdges = []; % edge IDs of each face                                                                                                                       
msh.faceElems = []; % left and right elements
msh.faceNums  = []; % left and right local face numbers
msh.faceR2L   = []; % face-to-face mapping index from right to left

% Augment metric information of faces for adaptive mesh
msh.faceLevel  = []; % tree level
msh.faceLID    = []; % ID among the leaf faces

% Augment topology information of faces for adaptive mesh
msh.faceParent   = []; % parent ID of each face
msh.faceChildren = []; % child IDs of each face

% Store other useful information for computation
msh.faceJac      = []; % Jacobian of the transformation from the referce face to a physical face
msh.faceIsCurved = []; % A face is indexed curved if any of its edge approximates a curved boundary
%**************************************************************************
% metric information of edges
msh.edgeType = []; % type

% topology information of edges
msh.edgeNodes  = []; % node IDs of each edge
msh.edgeNElems = []; % number of elements to which each edge belongs

% Augment topology information of edges for adaptive mesh
msh.edgeParent   = []; % parent ID of each edge
msh.edgeChildren = []; % child IDs of each edge

%**************************************************************************
% metric information of nodes
msh.nodeCoor = []; % coordinate of each node

end

