% Suppose the level difference between neighboring elements is no more than
% 1 when refining and coarsening
% Copyright (C) Xiangyi Meng.

function [msh, U] = refineAndCoarsenHexaMesh(msh, candReIDs, candCeIDs, nv, U, quad, bs, getFaceNode, getEdgeNode, layout)

if (nargin < 1)
    error('Not enough arguments')
end

if (msh.type ~= 312)
    error('The structure of the given mesh is not adaptive')
end

if (nargin < 2) 
    candReIDs = [];
end

if (nargin < 3) 
    candCeIDs = [];
end

if isempty(candReIDs) && isempty(candCeIDs)
    return
end

if (nargin < 4) || isempty(nv)
    nv = 1;
end

if (nargin < 5) || isempty(U)
    U = zeros(4, msh.nLElems);
end
if (size(U, 2) ~= msh.nLElems * nv)
    error('The given mesh and the solution coefficients does not coinside')
end

if (nargin < 6) || isempty(quad)
    quad = GaussQuadratureRule_hexa([2, 2, 2], 301);
end

if (nargin < 7) || isempty(bs)
    quad2 = GaussQuadratureRule_square([2, 2], 201);
    bs = setBasisFunctionSet_cube(quad2, quad, 301);
    bs = addBasisFunctionSet_cube(bs, quad2, quad);
end
if (size(U, 1) ~= bs.nb)
    error('The given U and bs does not coinside')
end

if (nargin < 8) || isempty(getFaceNode)
    getFaceNode = @(p1, p2, p3, p4, faceType)0.25 * (p1 + p2 + p3 + p4);
end

if (nargin < 9) || isempty(getEdgeNode)
    getEdgeNode = @(p1, p2, edgeType)0.5 * (p1 + p2);
end

if (nargin < 10) || isempty(layout)
    layout = 2;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

%**************************************************************************
% Firstly refine the elements. Note that we may add more elements to refine
% to ensure the level difference between neighboring elements is no more
% than 1
%**************************************************************************
% Delete the candidate elements whose level is greater or equal to maxLevel
reIDs = candReIDs(msh.elemLevel(candReIDs) < msh.maxLevel);
nre = length(reIDs);
shouldRefine = false(1, msh.nElems);
shouldRefine(reIDs) = true;

% Make sure the level difference of adjacent elements is no more than 1 
% by increasing the number of refined elements
reIDs(nre + 1 : floor(msh.nLElems / 3)) = 0; % preallocation
i = 1;
while i <= nre
    % fields of the current loop element
    ID = reIDs(i); 
    level     = msh.elemLevel(ID);
    faceIDs   = msh.elemFaces(:, ID);   
    faceTypes = msh.faceType(faceIDs);   
    faceElems = msh.faceElems(:, faceIDs);
    
    % Determine if refine the first adjacent element
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        ID1 = faceElems(1, 1) * (ID == faceElems(2, 1)) + faceElems(2, 1) * (ID == faceElems(1, 1));
        if (msh.elemLevel(ID1) < level) && ~shouldRefine(ID1)
            shouldRefine(ID1) = true;
            nre = nre + 1;
            reIDs(nre) = ID1;
        end        
    end
    % Determine if refine the second adjacent element
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        ID2 = faceElems(1, 2) * (ID == faceElems(2, 2)) + faceElems(2, 2) * (ID == faceElems(1, 2));
        if (msh.elemLevel(ID2) < level) && ~shouldRefine(ID2)
            shouldRefine(ID2) = true;
            nre = nre + 1;
            reIDs(nre) = ID2;
        end          
    end    
    % Determine if refine the third adjacent element
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        ID3 = faceElems(1, 3) * (ID == faceElems(2, 3)) + faceElems(2, 3) * (ID == faceElems(1, 3));
        if (msh.elemLevel(ID3) < level) && ~shouldRefine(ID3)
            shouldRefine(ID3) = true;
            nre = nre + 1;
            reIDs(nre) = ID3;
        end         
    end    
    % Determine if refine the fourth adjacent element
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        ID4 = faceElems(1, 4) * (ID == faceElems(2, 4)) + faceElems(2, 4) * (ID == faceElems(1, 4));
        if (msh.elemLevel(ID4) < level) && ~shouldRefine(ID4)
            shouldRefine(ID4) = true;
            nre = nre + 1;
            reIDs(nre) = ID4;
        end         
    end    
    % Determine if refine the fifth adjacent element
    if (faceTypes(5) == 0 || faceTypes(5) == 1)
        ID5 = faceElems(1, 5) * (ID == faceElems(2, 5)) + faceElems(2, 5) * (ID == faceElems(1, 5));
        if (msh.elemLevel(ID5) < level) && ~shouldRefine(ID5)
            shouldRefine(ID5) = true;
            nre = nre + 1;
            reIDs(nre) = ID5;
        end         
    end    
    % Determine if refine the sixth adjacent element
    if (faceTypes(6) == 0 || faceTypes(6) == 1)
        ID6 = faceElems(1, 6) * (ID == faceElems(2, 6)) + faceElems(2, 6) * (ID == faceElems(1, 6));
        if (msh.elemLevel(ID6) < level) && ~shouldRefine(ID6)
            shouldRefine(ID6) = true;
            nre = nre + 1;
            reIDs(nre) = ID6;
        end         
    end
    i = i + 1;
end
reIDs(nre + 1 : end) = [];

% We should first refine the elements of level 0, then level 1, then ...
% finally level (maxLvel - 1)
reLevels = msh.elemLevel(reIDs);
sreIDs = zeros(1, nre);
s = 0;
for i = 0 : msh.maxLevel - 1
    temp = reIDs(reLevels == i);
    l = length(temp);
    sreIDs(s + 1 : s + l) = temp;
    s = s + l;
end
if (s ~= nre)
    error('Wrong size after sorting')
end
reIDs = sreIDs;

% Then we refine the elements by adding new elements and faces
% preallocation for elements
msh.elemSize(msh.nElems + 1 : msh.nElems + 8 * nre)        = 0;
msh.elemDiam(msh.nElems + 1 : msh.nElems + 8 * nre)        = 0;
msh.elemNodes(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0;
msh.elemEdges(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0;
msh.elemFaces(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0;
msh.elemLevel(msh.nElems + 1 : msh.nElems + 8 * nre)       = 0;
msh.elemParent(msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemChildren(:, msh.nElems + 1 : msh.nElems + 8 * nre) = 0;
msh.elemGPx(:, msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemGPy(:, msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemGPz(:, msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemJac(:, msh.nElems + 1 : msh.nElems + 8 * nre)      = 0;
msh.elemJxix(:, msh.nElems + 1 : msh.nElems + 8 * nre)     = 0; 
msh.elemJxiy(:, msh.nElems + 1 : msh.nElems + 8 * nre)     = 0; 
msh.elemJxiz(:, msh.nElems + 1 : msh.nElems + 8 * nre)     = 0; 
msh.elemJetax(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0; 
msh.elemJetay(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0; 
msh.elemJetaz(:, msh.nElems + 1 : msh.nElems + 8 * nre)    = 0; 
msh.elemJzetax(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0; 
msh.elemJzetay(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0; 
msh.elemJzetaz(:, msh.nElems + 1 : msh.nElems + 8 * nre)   = 0;
msh.elemIsCurved(msh.nElems + 1 : msh.nElems + 8 * nre)    = false;
% preallocation for faces
msh.faceNormalx(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0; 
msh.faceNormaly(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0; 
msh.faceNormalz(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)  = 0;
msh.faceSize(msh.nFaces + 1 : msh.nFaces + 36 * nre)        = 0; 
msh.faceType(msh.nFaces + 1 : msh.nFaces + 36 * nre)        = 0;
msh.faceNodes(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)    = 0;
msh.faceEdges(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)    = 0;
msh.faceElems(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)    = 0; 
msh.faceNums(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)     = 0; 
msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 36 * nre)         = 0;
msh.faceLevel(msh.nFaces + 1 : msh.nFaces + 36 * nre)       = 0; 
msh.faceParent(msh.nFaces + 1 : msh.nFaces + 36 * nre)      = 0; 
msh.faceChildren(:, msh.nFaces + 1 : msh.nFaces + 36 * nre) = 0; 
msh.faceJac(:, msh.nFaces + 1 : msh.nFaces + 36 * nre)      = 0;  
msh.faceIsCurved(msh.nFaces + 1 : msh.nFaces + 36 * nre)    = false; 
% preallocation for edges
msh.edgeType(msh.nEdges + 1 : msh.nEdges + 66 * nre)        = 0; 
msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 66 * nre)    = 0; 
msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 66 * nre)      = 0; 
msh.edgeParent(msh.nEdges + 1 : msh.nEdges + 66 * nre)      = 0; 
msh.edgeChildren(:, msh.nEdges + 1 : msh.nEdges + 66 * nre) = 0;
% preallocation for nodes
msh.nodeCoor(:, msh.nNodes + 1 : msh.nNodes + 24 * nre) = 0; 
% preallocation for inverse of mass matrices
msh.massMatInv{msh.nElems + 1 + 8 * nre} = [];
% preallocation for solution
UT = zeros(size(U, 1), (msh.nElems + 8 * nre) * nv);
if (layout == 1)
    U = U(:, (0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
end
UT(:, (-nv + 1 : 0)' + nv * msh.LElems) = U;

ind = [1,2,4,5,10,11,13,14;2,3,5,6,11,12,14,15;5,6,8,9,14,15,17,18;4,5,7,8,13,14,16,17;10,11,13,14,19,20,22,23;11,12,14,15,20,21,23,24;14,15,17,18,23,24,26,27;13,14,16,17,22,23,25,26];
faceIDs1  = zeros(3, 2, 2);
faceIDs2  = zeros(2, 3, 2);
faceIDs3  = zeros(2, 2, 3);
edgeIDs1  = zeros(2, 3, 3);
edgeIDs2  = zeros(3, 2, 3);
edgeIDs3  = zeros(3, 3, 2);
nodeIDs   = zeros(3, 3, 3);
nodeCoors = zeros(3, 27);
for i = 1 : nre
    % current loop element and its nodes, faces and edges 
    ID = reIDs(i);
    elemNodes     = msh.elemNodes(:, ID);
    elemEdges     = msh.elemEdges(:, ID);
    elemFaces     = msh.elemFaces(:, ID);
    elemNodeCoors = msh.nodeCoor(:, elemNodes);
    elemEdgeTypes = msh.edgeType(:, elemEdges);
    elemEdgeNodes = msh.edgeNodes(:, elemEdges);
    
    nodeIDs([1, 3, 9, 7, 19, 21, 27, 25])      = elemNodes;
    nodeCoors(:, [1, 3, 9, 7, 19, 21, 27, 25]) = elemNodeCoors;
           
    % creation of subelements
    elemChildren = reshape(msh.nElems + 1 : msh.nElems + 8, [2, 2, 2]);
    msh.elemChildren(:, ID)      = elemChildren(:); 
    msh.elemLevel(elemChildren)  = msh.elemLevel(ID) + 1;
    msh.elemParent(elemChildren) = ID;   
    
    % creation of subedges and new nodes
    % 1st edge
    [nodeIDs(2), nodeCoors(:, 2), edgeChildren] = refineEdge(elemEdges(1), elemNodeCoors(:, [1, 2]), true);
    if (elemEdgeNodes(1, 1) == elemNodes(1))
        edgeIDs1(:, 1, 1) = edgeChildren;
    else
        edgeIDs1(:, 1, 1) = flip(edgeChildren);        
    end
    % 2nd edge
    [nodeIDs(8), nodeCoors(:, 8), edgeChildren] = refineEdge(elemEdges(2), elemNodeCoors(:, [4, 3]), true);
    if (elemEdgeNodes(1, 2) == elemNodes(4))
        edgeIDs1(:, 3, 1) = edgeChildren;
    else
        edgeIDs1(:, 3, 1) = flip(edgeChildren);        
    end       
    % 3rd edge
    [nodeIDs(20), nodeCoors(:, 20), edgeChildren] = refineEdge(elemEdges(3), elemNodeCoors(:, [5, 6]), true);
    if (elemEdgeNodes(1, 3) == elemNodes(5))
        edgeIDs1(:, 1, 3) = edgeChildren;
    else
        edgeIDs1(:, 1, 3) = flip(edgeChildren);        
    end      
    % 4th edge
    [nodeIDs(26), nodeCoors(:, 26), edgeChildren] = refineEdge(elemEdges(4), elemNodeCoors(:, [8, 7]), true);
    if (elemEdgeNodes(1, 4) == elemNodes(8))
        edgeIDs1(:, 3, 3) = edgeChildren;
    else
        edgeIDs1(:, 3, 3) = flip(edgeChildren);        
    end       
    % 5th edge
    [nodeIDs(4), nodeCoors(:, 4), edgeChildren] = refineEdge(elemEdges(5), elemNodeCoors(:, [1, 4]), true);
    if (elemEdgeNodes(1, 5) == elemNodes(1))
        edgeIDs2(1, :, 1) = edgeChildren;
    else
        edgeIDs2(1, :, 1) = flip(edgeChildren);        
    end     
    % 6th edge
    [nodeIDs(6), nodeCoors(:, 6), edgeChildren] = refineEdge(elemEdges(6), elemNodeCoors(:, [2, 3]), true);
    if (elemEdgeNodes(1, 6) == elemNodes(2))
        edgeIDs2(3, :, 1) = edgeChildren;
    else
        edgeIDs2(3, :, 1) = flip(edgeChildren);        
    end   
    % 7th edge
    [nodeIDs(22), nodeCoors(:, 22), edgeChildren] = refineEdge(elemEdges(7), elemNodeCoors(:, [5, 8]), true);
    if (elemEdgeNodes(1, 7) == elemNodes(5))
        edgeIDs2(1, :, 3) = edgeChildren;
    else
        edgeIDs2(1, :, 3) = flip(edgeChildren);        
    end  
    % 8th edge
    [nodeIDs(24), nodeCoors(:, 24), edgeChildren] = refineEdge(elemEdges(8), elemNodeCoors(:, [6, 7]), true);
    if (elemEdgeNodes(1, 8) == elemNodes(6))
        edgeIDs2(3, :, 3) = edgeChildren;
    else
        edgeIDs2(3, :, 3) = flip(edgeChildren);        
    end    
    % 9th edge
    [nodeIDs(10), nodeCoors(:, 10), edgeChildren] = refineEdge(elemEdges(9), elemNodeCoors(:, [1, 5]), true);
    if (elemEdgeNodes(1, 9) == elemNodes(1))
        edgeIDs3(1, 1, :) = edgeChildren;
    else
        edgeIDs3(1, 1, :) = flip(edgeChildren);        
    end 
    % 10th edge
    [nodeIDs(12), nodeCoors(:, 12), edgeChildren] = refineEdge(elemEdges(10), elemNodeCoors(:, [2, 6]), true);
    if (elemEdgeNodes(1, 10) == elemNodes(2))
        edgeIDs3(3, 1, :) = edgeChildren;
    else
        edgeIDs3(3, 1, :) = flip(edgeChildren);        
    end 
    % 11th edge
    [nodeIDs(16), nodeCoors(:, 16), edgeChildren] = refineEdge(elemEdges(11), elemNodeCoors(:, [4, 8]), true);
    if (elemEdgeNodes(1, 11) == elemNodes(4))
        edgeIDs3(1, 3, :) = edgeChildren;
    else
        edgeIDs3(1, 3, :) = flip(edgeChildren);        
    end    
    % 12th edge
    [nodeIDs(18), nodeCoors(:, 18), edgeChildren] = refineEdge(elemEdges(12), elemNodeCoors(:, [3, 7]), true);
    if (elemEdgeNodes(1, 12) == elemNodes(3))
        edgeIDs3(3, 3, :) = edgeChildren;
    else
        edgeIDs3(3, 3, :) = flip(edgeChildren);        
    end

    % creation of subfaces, new nodes and edges
    % 1st face
    faceEdgeIDs1 = reshape(edgeIDs3(1, :, :), [3, 2]);
    faceEdgeIDs2 = reshape(edgeIDs2(1, :, :), [2, 3]);
    faceNodeIDs  = reshape(nodeIDs(1, :, :), [1, 9]);
    faceElemIDs  = reshape(elemChildren(1, :, :), [1, 4]);
    [nodeIDs(13), nodeCoors(:, 13), edgeIDs, faceIDs] = refineFace(ID, elemFaces(1), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([9, 11, 5, 7]), elemNodeCoors(:, [1, 4, 8, 5]));
    edgeIDs3(1, 2, :) = edgeIDs(1 : 2);
    edgeIDs2(1, :, 2) = edgeIDs(3 : 4);
    faceIDs1(1, :, :) = reshape(faceIDs, [1, 2, 2]);
    % 2nd face
    faceEdgeIDs1 = reshape(edgeIDs3(3, :, :), [3, 2]);
    faceEdgeIDs2 = reshape(edgeIDs2(3, :, :), [2, 3]);
    faceNodeIDs  = reshape(nodeIDs(3, :, :), [1, 9]);
    faceElemIDs  = reshape(elemChildren(2, :, :), [1, 4]);
    [nodeIDs(15), nodeCoors(:, 15), edgeIDs, faceIDs] = refineFace(ID, elemFaces(2), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([10, 12, 6, 8]), elemNodeCoors(:, [2, 3, 7, 6]));
    edgeIDs3(3, 2, :) = edgeIDs(1 : 2);
    edgeIDs2(3, :, 2) = edgeIDs(3 : 4);
    faceIDs1(3, :, :) = reshape(faceIDs, [1, 2, 2]);    
    % 3rd face
    faceEdgeIDs1 = reshape(edgeIDs1([1, 7, 13, 2, 8, 14]), [3, 2]);
    faceEdgeIDs2 = reshape(edgeIDs3([1, 10, 2, 11, 3, 12]), [2, 3]);
    faceNodeIDs  = nodeIDs([1, 10, 19, 2, 11, 20, 3, 12, 21]);
    faceElemIDs  = elemChildren([1, 5, 2, 6]);
    [nodeIDs(11), nodeCoors(:, 11), edgeIDs, faceIDs] = refineFace(ID, elemFaces(3), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([1, 3, 9, 10]), elemNodeCoors(:, [1, 5, 6, 2]));
    edgeIDs1(:, 1, 2) = edgeIDs(1 : 2);
    edgeIDs3(2, 1, :) = edgeIDs(3 : 4);
    faceIDs2(:, 1, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);
    % 4th face
    faceEdgeIDs1 = reshape(edgeIDs1([5, 11, 17, 6, 12, 18]), [3, 2]);
    faceEdgeIDs2 = reshape(edgeIDs3([7, 16, 8, 17, 9, 18]), [2, 3]);
    faceNodeIDs  = nodeIDs([7, 16, 25, 8, 17, 26, 9, 18, 27]);
    faceElemIDs  = elemChildren([3, 7, 4, 8]);
    [nodeIDs(17), nodeCoors(:, 17), edgeIDs, faceIDs] = refineFace(ID, elemFaces(4), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([2, 4, 11, 12]), elemNodeCoors(:, [4, 8, 7, 3]));
    edgeIDs1(:, 3, 2) = edgeIDs(1 : 2);
    edgeIDs3(2, 3, :) = edgeIDs(3 : 4);
    faceIDs2(:, 3, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);  
    % 5th face
    faceEdgeIDs1 = edgeIDs2(:, :, 1);
    faceEdgeIDs2 = edgeIDs1(:, :, 1);
    faceNodeIDs  = reshape(nodeIDs(:, :, 1), [1, 9]);
    faceElemIDs  = reshape(elemChildren(:, :, 1), [1, 4]);
    [nodeIDs(5), nodeCoors(:, 5), edgeIDs, faceIDs] = refineFace(ID, elemFaces(5), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([5, 6, 1, 2]), elemNodeCoors(:, 1 : 4));
    edgeIDs2(2, :, 1) = edgeIDs(1 : 2);
    edgeIDs1(:, 2, 1) = edgeIDs(3 : 4);
    faceIDs3(:, :, 1) = reshape(faceIDs, [2, 2]);   
    % 6th face
    faceEdgeIDs1 = edgeIDs2(:, :, 3);
    faceEdgeIDs2 = edgeIDs1(:, :, 3);
    faceNodeIDs  = reshape(nodeIDs(:, :, 3), [1, 9]);
    faceElemIDs  = reshape(elemChildren(:, :, 2), [1, 4]);
    [nodeIDs(23), nodeCoors(:, 23), edgeIDs, faceIDs] = refineFace(ID, elemFaces(6), faceEdgeIDs1, faceEdgeIDs2, faceNodeIDs, faceElemIDs, elemEdgeTypes([7, 8, 3, 4]), elemNodeCoors(:, 5 : 8));
    edgeIDs2(2, :, 3) = edgeIDs(1 : 2);
    edgeIDs1(:, 2, 3) = edgeIDs(3 : 4);
    faceIDs3(:, :, 3) = reshape(faceIDs, [2, 2]);  
    % interior nodes, edges and faces
    nodeIDs(14) = msh.nNodes + 1;
    nodeCoors(:, 14) = sum(nodeCoors(:, [5, 11, 13, 15, 17, 23]), 2) / 6;
    msh.nodeCoor(:, nodeIDs(14)) = nodeCoors(:, 14);
    msh.nNodes = msh.nNodes + 1;
    
    edgeIDs1(:, 2, 2) = msh.nEdges + 1 : msh.nEdges + 2;
    msh.edgeType(msh.nEdges + 1 : msh.nEdges + 2)        = 0;
    msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2)    = [nodeIDs(13), nodeIDs(14); nodeIDs(14), nodeIDs(15)];
    msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)      = 4;  
    msh.edgeParent(msh.nEdges + 1 : msh.nEdges + 2)      = 0;
    msh.edgeChildren(:, msh.nEdges + 1 : msh.nEdges + 2) = 0;
    msh.nEdges = msh.nEdges + 2;   
    
    edgeIDs2(2, :, 2) = msh.nEdges + 1 : msh.nEdges + 2;
    msh.edgeType(msh.nEdges + 1 : msh.nEdges + 2)        = 0;
    msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2)    = [nodeIDs(11), nodeIDs(14); nodeIDs(14), nodeIDs(17)];
    msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)      = 4;  
    msh.edgeParent(msh.nEdges + 1 : msh.nEdges + 2)      = 0;
    msh.edgeChildren(:, msh.nEdges + 1 : msh.nEdges + 2) = 0;
    msh.nEdges = msh.nEdges + 2;   
    
    edgeIDs3(2, 2, :) = msh.nEdges + 1 : msh.nEdges + 2;
    msh.edgeType(msh.nEdges + 1 : msh.nEdges + 2)        = 0;
    msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2)    = [nodeIDs(5), nodeIDs(14); nodeIDs(14), nodeIDs(23)];
    msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)      = 4; 
    msh.edgeParent(msh.nEdges + 1 : msh.nEdges + 2)      = 0;
    msh.edgeChildren(:, msh.nEdges + 1 : msh.nEdges + 2) = 0;
    msh.nEdges = msh.nEdges + 2;     
    
    faceIDs1(2, :, :) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [1, 2, 2]);
    msh.faceNodes(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 1 : 2, 1 : 2), [1, 4]);
    msh.faceNodes(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 2 : 3, 1 : 2), [1, 4]);
    msh.faceNodes(3, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 2 : 3, 2 : 3), [1, 4]);
    msh.faceNodes(4, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 1 : 2, 2 : 3), [1, 4]);  
    msh.faceEdges(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs3(2, 1 : 2, :), [1, 4]);
    msh.faceEdges(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs3(2, 2 : 3, :), [1, 4]);
    msh.faceEdges(3, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs2(2, :, 1 : 2), [1, 4]);
    msh.faceEdges(4, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs2(2, :, 2 : 3), [1, 4]);    
    msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemChildren(1, :, :), [1, 4]);
    msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemChildren(2, :, :), [1, 4]);
    msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 2;
    msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 1;
    msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4)      = 1;
    tempNormals = 0.5 * mycross3(nodeCoors(:, [2, 5, 11, 14]) - nodeCoors(:, [14, 17, 23, 26]), nodeCoors(:, [5, 8, 14, 17]) - nodeCoors(:, [11, 14, 20, 23]));
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    msh.faceNormalx(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(1, :);
    msh.faceNormaly(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(2, :);
    msh.faceNormalz(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(3, :);
    msh.faceSize(msh.nFaces + 1 : msh.nFaces + 4)    = tempSizes;
    msh.faceJac(msh.nFaces + 1 : msh.nFaces + 4)     = tempSizes / 4;
    msh.nFaces = msh.nFaces + 4; 
    
    faceIDs2(:, 2, :) = reshape(msh.nFaces + [1, 3, 2, 4], [2, 1, 2]); 
    msh.faceNodes(1, msh.nFaces + 1 : msh.nFaces + 4) = nodeIDs([4, 13, 5, 14]);
    msh.faceNodes(2, msh.nFaces + 1 : msh.nFaces + 4) = nodeIDs([13, 22, 14, 23]);
    msh.faceNodes(3, msh.nFaces + 1 : msh.nFaces + 4) = nodeIDs([14, 23, 15, 24]);
    msh.faceNodes(4, msh.nFaces + 1 : msh.nFaces + 4) = nodeIDs([5, 14, 6, 15]);  
    msh.faceEdges(1, msh.nFaces + 1 : msh.nFaces + 4) = edgeIDs1([3, 9, 4, 10]);
    msh.faceEdges(2, msh.nFaces + 1 : msh.nFaces + 4) = edgeIDs1([9, 15, 10, 16]);
    msh.faceEdges(3, msh.nFaces + 1 : msh.nFaces + 4) = edgeIDs3([4, 13, 5, 14]);
    msh.faceEdges(4, msh.nFaces + 1 : msh.nFaces + 4) = edgeIDs3([5, 14, 6, 15]);    
    msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = elemChildren([1, 5, 2, 6]);
    msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = elemChildren([3, 7, 4, 8]);
    msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 4;
    msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 3;    
    msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4) = 1;
    tempNormals = 0.5 * mycross3(nodeCoors(:, [4, 13, 5, 14]) - nodeCoors(:, [14, 23, 15, 24]), nodeCoors(:, [13, 22, 14, 23]) - nodeCoors(:, [5, 14, 6, 15]));
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;   
    msh.faceNormalx(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(1, :);
    msh.faceNormaly(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(2, :);
    msh.faceNormalz(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(3, :);
    msh.faceSize(msh.nFaces + 1 : msh.nFaces + 4)    = tempSizes;
    msh.faceJac(msh.nFaces + 1 : msh.nFaces + 4)     = tempSizes / 4;    
    msh.nFaces = msh.nFaces + 4;
    
    faceIDs3(:, :, 2) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [2, 2]);
    msh.faceNodes(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(1 : 2, 1 : 2, 2), [1, 4]);
    msh.faceNodes(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2 : 3, 1 : 2, 2), [1, 4]);
    msh.faceNodes(3, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2 : 3, 2 : 3, 2), [1, 4]);
    msh.faceNodes(4, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(1 : 2, 2 : 3, 2), [1, 4]); 
    msh.faceEdges(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs2(1 : 2, :, 2), [1, 4]);
    msh.faceEdges(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs2(2 : 3, :, 2), [1, 4]);
    msh.faceEdges(3, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs1(:, 1 : 2, 2), [1, 4]);
    msh.faceEdges(4, msh.nFaces + 1 : msh.nFaces + 4) = reshape(edgeIDs1(:, 2 : 3, 2), [1, 4]);    
    msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemChildren(:, :, 1), [1, 4]);
    msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemChildren(:, :, 2), [1, 4]);
    msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 6;
    msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 5;    
    msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4) = 1;
    tempNormals = 0.5 * mycross3(nodeCoors(:, [10, 11, 13, 14]) - nodeCoors(:, [14, 15, 17, 18]), nodeCoors(:, [11, 12, 14, 15]) - nodeCoors(:, [13, 14, 16, 17]));
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes; 
    msh.faceNormalx(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(1, :);
    msh.faceNormaly(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(2, :);
    msh.faceNormalz(msh.nFaces + 1 : msh.nFaces + 4) = tempNormals(3, :);
    msh.faceSize(msh.nFaces + 1 : msh.nFaces + 4)    = tempSizes;
    msh.faceJac(msh.nFaces + 1 : msh.nFaces + 4)     = tempSizes / 4;    
    msh.nFaces = msh.nFaces + 4;
    
    % element information 
    subElemNodeCoorx = reshape(nodeCoors(1, ind), [8, 8]);
    subElemNodeCoory = reshape(nodeCoors(2, ind), [8, 8]);
    subElemNodeCoorz = reshape(nodeCoors(3, ind), [8, 8]);
    
    x     = quad.N      * subElemNodeCoorx;
    xxi   = quad.N_xi   * subElemNodeCoorx;
    xeta  = quad.N_eta  * subElemNodeCoorx;
    xzeta = quad.N_zeta * subElemNodeCoorx;
    y     = quad.N      * subElemNodeCoory;
    yxi   = quad.N_xi   * subElemNodeCoory;
    yeta  = quad.N_eta  * subElemNodeCoory;   
    yzeta = quad.N_zeta * subElemNodeCoory;  
    z     = quad.N      * subElemNodeCoorz;
    zxi   = quad.N_xi   * subElemNodeCoorz;
    zeta  = quad.N_eta  * subElemNodeCoorz;
    zzeta = quad.N_zeta * subElemNodeCoorz;  
    J     = xxi .* (yeta .* zzeta - yzeta .* zeta) + xeta .* (yzeta .* zxi - yxi .* zzeta) + xzeta .* (yxi .* zeta - yeta .* zxi);
    
    si = 1 / 6 * dot(nodeCoors(:, ind(7, :)) - nodeCoors(:, ind(1, :)), mycross3(nodeCoors(:, ind(2, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(3, :)) - nodeCoors(:, ind(6, :))) + mycross3(nodeCoors(:, ind(4, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(8, :)) - nodeCoors(:, ind(3, :))) + mycross3(nodeCoors(:, ind(5, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(6, :)) - nodeCoors(:, ind(8, :))));
    msh.elemSize(elemChildren)      = si;
    msh.elemDiam(elemChildren)      = si.^(1/3);
    msh.elemNodes(1, elemChildren)  = reshape(nodeIDs(1 : 2, 1 : 2, 1 : 2), [1, 8]);
    msh.elemNodes(2, elemChildren)  = reshape(nodeIDs(2 : 3, 1 : 2, 1 : 2), [1, 8]);
    msh.elemNodes(3, elemChildren)  = reshape(nodeIDs(2 : 3, 2 : 3, 1 : 2), [1, 8]);
    msh.elemNodes(4, elemChildren)  = reshape(nodeIDs(1 : 2, 2 : 3, 1 : 2), [1, 8]);
    msh.elemNodes(5, elemChildren)  = reshape(nodeIDs(1 : 2, 1 : 2, 2 : 3), [1, 8]);
    msh.elemNodes(6, elemChildren)  = reshape(nodeIDs(2 : 3, 1 : 2, 2 : 3), [1, 8]);
    msh.elemNodes(7, elemChildren)  = reshape(nodeIDs(2 : 3, 2 : 3, 2 : 3), [1, 8]);
    msh.elemNodes(8, elemChildren)  = reshape(nodeIDs(1 : 2, 2 : 3, 2 : 3), [1, 8]);
    msh.elemEdges(1, elemChildren)  = reshape(edgeIDs1(:, 1 : 2, 1 : 2), [1, 8]);
    msh.elemEdges(2, elemChildren)  = reshape(edgeIDs1(:, 2 : 3, 1 : 2), [1, 8]);
    msh.elemEdges(3, elemChildren)  = reshape(edgeIDs1(:, 1 : 2, 2 : 3), [1, 8]);
    msh.elemEdges(4, elemChildren)  = reshape(edgeIDs1(:, 2 : 3, 2 : 3), [1, 8]);  
    msh.elemEdges(5, elemChildren)  = reshape(edgeIDs2(1 : 2, :, 1 : 2), [1, 8]);
    msh.elemEdges(6, elemChildren)  = reshape(edgeIDs2(2 : 3, :, 1 : 2), [1, 8]);
    msh.elemEdges(7, elemChildren)  = reshape(edgeIDs2(1 : 2, :, 2 : 3), [1, 8]);
    msh.elemEdges(8, elemChildren)  = reshape(edgeIDs2(2 : 3, :, 2 : 3), [1, 8]);    
    msh.elemEdges(9, elemChildren)  = reshape(edgeIDs3(1 : 2, 1 : 2, :), [1, 8]);
    msh.elemEdges(10, elemChildren) = reshape(edgeIDs3(2 : 3, 1 : 2, :), [1, 8]);
    msh.elemEdges(11, elemChildren) = reshape(edgeIDs3(1 : 2, 2 : 3, :), [1, 8]);
    msh.elemEdges(12, elemChildren) = reshape(edgeIDs3(2 : 3, 2 : 3, :), [1, 8]);    
    msh.elemFaces(1, elemChildren)  = reshape(faceIDs1(1 : 2, :, :), [1, 8]);
    msh.elemFaces(2, elemChildren)  = reshape(faceIDs1(2 : 3, :, :), [1, 8]);
    msh.elemFaces(3, elemChildren)  = reshape(faceIDs2(:, 1 : 2, :), [1, 8]);
    msh.elemFaces(4, elemChildren)  = reshape(faceIDs2(:, 2 : 3, :), [1, 8]);
    msh.elemFaces(5, elemChildren)  = reshape(faceIDs3(:, :, 1 : 2), [1, 8]);
    msh.elemFaces(6, elemChildren)  = reshape(faceIDs3(:, :, 2 : 3), [1, 8]); 
    msh.elemGPx(:, elemChildren)    = x;
    msh.elemGPy(:, elemChildren)    = y;
    msh.elemGPz(:, elemChildren)    = z;
    msh.elemJac(:, elemChildren)    = J;
    msh.elemJxix(:, elemChildren)   = yeta .* zzeta - yzeta .* zeta;
    msh.elemJxiy(:, elemChildren)   = zeta .* xzeta - zzeta .* xeta;
    msh.elemJxiz(:, elemChildren)   = xeta .* yzeta - xzeta .* yeta;
    msh.elemJetax(:, elemChildren)  = yzeta .* zxi - yxi .* zzeta;
    msh.elemJetay(:, elemChildren)  = zzeta .* xxi - zxi .* xzeta;
    msh.elemJetaz(:, elemChildren)  = xzeta .* yxi - xxi .* yzeta;
    msh.elemJzetax(:, elemChildren) = yxi .* zeta - yeta .* zxi;
    msh.elemJzetay(:, elemChildren) = zxi .* xeta - zeta .* xxi;
    msh.elemJzetaz(:, elemChildren) = xxi .* yeta - xeta .* yxi;
    msh.nElems = msh.nElems + 8;      
    
    % inverse of mass matrix
    mi1 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 1)));
    mi2 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 2)));
    mi3 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 3)));
    mi4 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 4)));
    mi5 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 5)));
    mi6 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 6)));
    mi7 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 7)));
    mi8 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 8)));
    msh.massMatInv{elemChildren(1)} = mi1;
    msh.massMatInv{elemChildren(2)} = mi2;
    msh.massMatInv{elemChildren(3)} = mi3;
    msh.massMatInv{elemChildren(4)} = mi4;
    msh.massMatInv{elemChildren(5)} = mi5;
    msh.massMatInv{elemChildren(6)} = mi6;
    msh.massMatInv{elemChildren(7)} = mi7;
    msh.massMatInv{elemChildren(8)} = mi8;
    
    % Project the solution in the coarse element to the eight refined
    % elements
    UE = UT(:, nv * (ID - 1) + 1 : nv * ID);
    if msh.elemIsCurved(ID)
        waring('Not implemented for curved elements yet')
        UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) = mi1 * (bs.phitw{1} * (bs.phi_sube{1} * UE .* J(:, 1)));   
        UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) = mi2 * (bs.phitw{1} * (bs.phi_sube{2} * UE .* J(:, 2)));
        UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) = mi3 * (bs.phitw{1} * (bs.phi_sube{3} * UE .* J(:, 3)));
        UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4)) = mi4 * (bs.phitw{1} * (bs.phi_sube{4} * UE .* J(:, 4))); 
        UT(:, nv * (elemChildren(5) - 1) + 1 : nv * elemChildren(5)) = mi5 * (bs.phitw{1} * (bs.phi_sube{5} * UE .* J(:, 5)));   
        UT(:, nv * (elemChildren(6) - 1) + 1 : nv * elemChildren(6)) = mi6 * (bs.phitw{1} * (bs.phi_sube{6} * UE .* J(:, 6)));
        UT(:, nv * (elemChildren(7) - 1) + 1 : nv * elemChildren(7)) = mi7 * (bs.phitw{1} * (bs.phi_sube{7} * UE .* J(:, 7)));
        UT(:, nv * (elemChildren(8) - 1) + 1 : nv * elemChildren(8)) = mi8 * (bs.phitw{1} * (bs.phi_sube{8} * UE .* J(:, 8)));          
    else        
        UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) = mi1 * (bs.phitw{1} * (bs.phi_sube{1} * UE .* J(:, 1)));   
        UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) = mi2 * (bs.phitw{1} * (bs.phi_sube{2} * UE .* J(:, 2)));
        UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) = mi3 * (bs.phitw{1} * (bs.phi_sube{3} * UE .* J(:, 3)));
        UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4)) = mi4 * (bs.phitw{1} * (bs.phi_sube{4} * UE .* J(:, 4))); 
        UT(:, nv * (elemChildren(5) - 1) + 1 : nv * elemChildren(5)) = mi5 * (bs.phitw{1} * (bs.phi_sube{5} * UE .* J(:, 5)));   
        UT(:, nv * (elemChildren(6) - 1) + 1 : nv * elemChildren(6)) = mi6 * (bs.phitw{1} * (bs.phi_sube{6} * UE .* J(:, 6)));
        UT(:, nv * (elemChildren(7) - 1) + 1 : nv * elemChildren(7)) = mi7 * (bs.phitw{1} * (bs.phi_sube{7} * UE .* J(:, 7)));
        UT(:, nv * (elemChildren(8) - 1) + 1 : nv * elemChildren(8)) = mi8 * (bs.phitw{1} * (bs.phi_sube{8} * UE .* J(:, 8)));   
    end
end
msh.faceNormalx(:, msh.nFaces + 1 : end)  = []; 
msh.faceNormaly(:, msh.nFaces + 1 : end)  = [];
msh.faceNormalz(:, msh.nFaces + 1 : end)  = []; 
msh.faceSize(msh.nFaces + 1 : end)        = []; 
msh.faceType(msh.nFaces + 1 : end)        = [];
msh.faceNodes(:, msh.nFaces + 1 : end)    = [];
msh.faceEdges(:, msh.nFaces + 1 : end)    = [];
msh.faceElems(:, msh.nFaces + 1 : end)    = []; 
msh.faceNums(:, msh.nFaces + 1 : end)     = []; 
msh.faceR2L(msh.nFaces + 1 : end)         = [];
msh.faceLevel(msh.nFaces + 1 : end)       = []; 
msh.faceParent(msh.nFaces + 1 : end)      = []; 
msh.faceChildren(:, msh.nFaces + 1 : end) = []; 
msh.faceJac(:, msh.nFaces + 1 : end)      = []; 
msh.faceIsCurved(msh.nFaces + 1 : end)    = []; 
msh.edgeType(msh.nEdges + 1 : end)        = []; 
msh.edgeNodes(:, msh.nEdges + 1 : end)    = []; 
msh.edgeNElems(msh.nEdges + 1 : end)      = []; 
msh.edgeParent(msh.nEdges + 1 : end)      = []; 
msh.edgeChildren(:, msh.nEdges + 1 : end) = [];
msh.nodeCoor(:, msh.nNodes + 1 : end)     = [];
msh.massMatInv(msh.nElems + 1 : end)      = [];

%**************************************************************************
% Then coarsen the elements. Note that we may delete some elements to
% coarsen to ensure the level difference between neighboring elements is no 
% more than 1
%**************************************************************************
% Delete the candidate elements which don't have a parent or should refine
ceIDs = candCeIDs(msh.elemParent(candCeIDs) ~= 0 & ~shouldRefine(candCeIDs));
ceIDs = sort(ceIDs);
nce = length(ceIDs);
shouldCoarsen = false(1, msh.nElems);
shouldCoarsen(ceIDs) = true;

% Make sure all the eight siblings need to be coarsened and store the IDs of
% their parents
npe = 0;
peIDs = zeros(1, floor(nce / 8));
i = 1;
while i <= nce
    ID = ceIDs(i);
    parentID = msh.elemParent(ID);
    % Make sure all the eight siblings need to be coarsened
    siblingIDs = msh.elemChildren(:, parentID);
    if any(~shouldCoarsen(siblingIDs))
        shouldCoarsen(ID) = false;
        i = i + 1;
        continue;
    end
    % Store the parent ID
    npe = npe + 1;
    peIDs(npe) = parentID;
    
    i = i + 8;
end
peIDs(npe + 1 : end) = [];

% We should first coarsen the elements of level maxLevel, then level 
% maxLevel - 1, then ..., finally level 1
peLevels = msh.elemLevel(peIDs);
speIDs = zeros(1, npe);
s = 0;
for i = msh.maxLevel - 1 : -1 : 0
    temp = peIDs(peLevels == i);
    l = length(temp);
    speIDs(s + 1 : s + l) = temp;
    s = s + l;
end
if (s ~= npe)
    error('Wrong size after sorting')
end
peIDs = speIDs;

% Update face-element connections, or delete elements, faces, edges, nodes 
% and store their IDs 
deIDs  = zeros(1,  8 * npe);
dfIDs  = zeros(1, 36 * npe); 
dedIDs = zeros(1, 66 * npe); 
dnIDs  = zeros(1, 24 * npe); 
nde  = 0;
ndf  = 0;
nded = 0;
ndn  = 0;
for i = 1 : npe
    ID = peIDs(i);
    elemNodes      = msh.elemNodes(:, ID);
    elemEdges      = msh.elemEdges(:, ID);
    elemFaces      = msh.elemFaces(:, ID);
    elemChildren   = msh.elemChildren(:, ID);
    elemChildEdges = msh.elemEdges(:, elemChildren);
    elemChildNodes = msh.elemNodes(:, elemChildren);
    
    % Make sure the level difference of adjacent elements is less than 1
    if any(msh.faceChildren(1, msh.faceChildren(:, elemFaces)) ~= 0)
        shouldCoarsen(elemChildren) = false;
        continue;  
    end
    
    % Project the solution in the eight refined elements into the coarse
    % element
    J = msh.elemJac(:, elemChildren);
    u8 = bs.phi{1} * UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(8));
    if msh.elemIsCurved(ID) 
        waring('Not implemented for curved elements yet')
        UT(:, nv * (ID - 1) + 1 : nv * ID) = msh.massMatInv{ID} * (bs.phitw_sube{1} * (u8(:, 1 : nv) .* J(:, 1)) ...
                                                                 + bs.phitw_sube{2} * (u8(:, nv + 1 : 2 * nv) .* J(:, 2)) ...
                                                                 + bs.phitw_sube{3} * (u8(:, 2 * nv + 1 : 3 * nv) .* J(:, 3)) ...
                                                                 + bs.phitw_sube{4} * (u8(:, 3 * nv + 1 : 4 * nv) .* J(:, 4)) ...
                                                                 + bs.phitw_sube{5} * (u8(:, 4 * nv + 1 : 5 * nv) .* J(:, 5)) ...
                                                                 + bs.phitw_sube{6} * (u8(:, 5 * nv + 1 : 6 * nv) .* J(:, 6)) ...
                                                                 + bs.phitw_sube{7} * (u8(:, 6 * nv + 1 : 7 * nv) .* J(:, 7)) ...
                                                                 + bs.phitw_sube{8} * (u8(:, 7 * nv + 1 : end) .* J(:, 8)));        
    else
        UT(:, nv * (ID - 1) + 1 : nv * ID) = msh.massMatInv{ID} * (bs.phitw_sube{1} * (u8(:, 1 : nv) .* J(:, 1)) ...
                                                                 + bs.phitw_sube{2} * (u8(:, nv + 1 : 2 * nv) .* J(:, 2)) ...
                                                                 + bs.phitw_sube{3} * (u8(:, 2 * nv + 1 : 3 * nv) .* J(:, 3)) ...
                                                                 + bs.phitw_sube{4} * (u8(:, 3 * nv + 1 : 4 * nv) .* J(:, 4)) ...
                                                                 + bs.phitw_sube{5} * (u8(:, 4 * nv + 1 : 5 * nv) .* J(:, 5)) ...
                                                                 + bs.phitw_sube{6} * (u8(:, 5 * nv + 1 : 6 * nv) .* J(:, 6)) ...
                                                                 + bs.phitw_sube{7} * (u8(:, 6 * nv + 1 : 7 * nv) .* J(:, 7)) ...
                                                                 + bs.phitw_sube{8} * (u8(:, 7 * nv + 1 : end) .* J(:, 8))); 
    end
                                                 
    % Update face-element connections or delete subfaces and associated 
    % edges and nodes 
    % 1st face
    edgeIDs = [elemChildEdges(11, [1, 5]), elemChildEdges(7, [1, 3])];
    nodeIDs = [elemNodes([1, 4, 8, 5])', elemChildNodes(8, 1)];
    coarsenFace(ID, elemFaces(1), edgeIDs, nodeIDs);
    % 2nd face
    edgeIDs = [elemChildEdges(12, [2, 6]), elemChildEdges(8, [2, 4])];
    nodeIDs = [elemNodes([2, 3, 7, 6])', elemChildNodes(2, 8)];
    coarsenFace(ID, elemFaces(2), edgeIDs, nodeIDs);
    % 3rd face
    edgeIDs = [elemChildEdges(3, [1, 2]), elemChildEdges(10, [1, 5])];
    nodeIDs = [elemNodes([1, 5, 6, 2])', elemChildNodes(6, 1)];
    coarsenFace(ID, elemFaces(3), edgeIDs, nodeIDs);
    % 4th face
    edgeIDs = [elemChildEdges(4, [3, 4]), elemChildEdges(12, [3, 7])];
    nodeIDs = [elemNodes([4, 8, 7, 3])', elemChildNodes(4, 8)];
    coarsenFace(ID, elemFaces(4), edgeIDs, nodeIDs);
    % 5th face
    edgeIDs = [elemChildEdges(6, [1, 3]), elemChildEdges(2, [1, 2])];
    nodeIDs = [elemNodes(1 : 4)', elemChildNodes(3, 1)];
    coarsenFace(ID, elemFaces(5), edgeIDs, nodeIDs);
    % 6th face
    edgeIDs = [elemChildEdges(8, [5, 7]), elemChildEdges(4, [5, 6])];
    nodeIDs = [elemNodes(5 : 8)', elemChildNodes(5, 8)];
    coarsenFace(ID, elemFaces(6), edgeIDs, nodeIDs);
    
    % Deal with the edges and nodes on the boundary
    subedgeIDs = msh.edgeChildren(:, elemEdges);
    nae = msh.edgeNElems(subedgeIDs) - 1;
    msh.edgeNElems(subedgeIDs) = nae;
    
    delSubedgeIDs = subedgeIDs(:, nae(1 : 2 : end) == 0);
    if ~isempty(delSubedgeIDs)
        ndse = numel(delSubedgeIDs);
        dedIDs(nded + 1 : nded + ndse) = delSubedgeIDs;
        nded = nded + ndse;
        
        dnIDs(ndn + 1 : ndn + 0.5 * ndse) = msh.edgeNodes(2, delSubedgeIDs(1, :));
        ndn = ndn + 0.5 * ndse;
    end
    
    % Delete interior faces, edges and nodes
    elemChildFace = msh.elemFaces(2, elemChildren(1));
    dfIDs(ndf + 1 : ndf + 12) = elemChildFace : elemChildFace + 11;
    ndf = ndf + 12;
    
    dedIDs(nded + 1 : nded + 6) = elemChildEdges(4, 1) : elemChildEdges(4, 1) + 5;
    nded = nded + 6;
    
    dnIDs(ndn + 1) = elemChildNodes(7, 1);
    ndn = ndn + 1;
    
    % delete elements
    msh.elemChildren(:, ID) = 0;
    deIDs(nde + 1 : nde + 8) = elemChildren;
    nde = nde + 8;    
end
deIDs(nde + 1 : end)   = [];
dfIDs(ndf + 1 : end)   = [];
dedIDs(nded + 1 : end) = [];
dnIDs(ndn + 1 : end)   = [];

% Re-index the remaining elements, faces, edges and nodes
remEIDs  = 1 : msh.nElems;
remFIDs  = 1 : msh.nFaces;
remEDIDs = 1 : msh.nEdges;
remNIDs  = 1 : msh.nNodes;
remEIDs(deIDs)   = []; % remaining element old IDs
remFIDs(dfIDs)   = []; % remaining face old IDs
remEDIDs(dedIDs) = []; % remaining edge old IDs
remNIDs(dnIDs)   = []; % remaining node old IDs

neIDs  = zeros(1, msh.nElems); 
nfIDs  = zeros(1, msh.nFaces); 
nedIDs = zeros(1, msh.nEdges);
nnIDs  = zeros(1, msh.nNodes);
neIDs(remEIDs)   = 1 : length(remEIDs);  % new element IDs
nfIDs(remFIDs)   = 1 : length(remFIDs);  % new face IDs
nedIDs(remEDIDs) = 1 : length(remEDIDs); % new element IDs
nnIDs(remNIDs)   = 1 : length(remNIDs);  % new node IDs

% Refresh element-node, element-edge, element-face connections, parent and 
% child IDs of elements
msh.elemNodes(:, remEIDs) = nnIDs(msh.elemNodes(:, remEIDs));
msh.elemEdges(:, remEIDs) = nedIDs(msh.elemEdges(:, remEIDs));
msh.elemFaces(:, remEIDs) = nfIDs(msh.elemFaces(:, remEIDs));

remEIDs_hasParent = remEIDs(msh.elemParent(remEIDs) ~= 0);
msh.elemParent(remEIDs_hasParent) = neIDs(msh.elemParent(remEIDs_hasParent));

remEIDs_hasChild = remEIDs(msh.elemChildren(1, remEIDs) ~= 0);
msh.elemChildren(:, remEIDs_hasChild) = neIDs(msh.elemChildren(:, remEIDs_hasChild));

% Refresh face-element, face-node, face-edge connections, parent and child 
% IDs of faces
msh.faceNodes(:, remFIDs) = nnIDs(msh.faceNodes(:, remFIDs));

msh.faceEdges(:, remFIDs) = nedIDs(msh.faceEdges(:, remFIDs));

msh.faceElems(1, remFIDs) = neIDs(msh.faceElems(1, remFIDs));
remFIDs_int = remFIDs(msh.faceType(remFIDs) == 0 | msh.faceType(remFIDs) == 1);
msh.faceElems(2, remFIDs_int) = neIDs(msh.faceElems(2, remFIDs_int));

remFIDs_hasParent = remFIDs(msh.faceParent(remFIDs) ~= 0);
msh.faceParent(remFIDs_hasParent) = nfIDs(msh.faceParent(remFIDs_hasParent));

remFIDs_hasChild = remFIDs(msh.faceChildren(1, remFIDs) ~= 0);
msh.faceChildren(:, remFIDs_hasChild) = nfIDs(msh.faceChildren(:, remFIDs_hasChild));

% Refresh edge-node connections, parent and child IDs of edges
msh.edgeNodes(:, remEDIDs) = nnIDs(msh.edgeNodes(:, remEDIDs));

remEDIDs_hasParent = remEDIDs(msh.edgeParent(remEDIDs) ~= 0);
msh.edgeParent(remEDIDs_hasParent) = nedIDs(msh.edgeParent(remEDIDs_hasParent));

remEDIDs_hasChild = remEDIDs(msh.edgeChildren(1, remEDIDs) ~= 0);
msh.edgeChildren(:, remEDIDs_hasChild) = nedIDs(msh.edgeChildren(:, remEDIDs_hasChild));

% Finally remove the deleted elements and faces from the mesh
msh.elemSize(deIDs)        = [];
msh.elemDiam(deIDs)        = [];
msh.elemNodes(:, deIDs)    = [];
msh.elemEdges(:, deIDs)    = [];
msh.elemFaces(:, deIDs)    = [];
msh.elemLevel(deIDs)       = [];
msh.elemParent(deIDs)      = [];
msh.elemChildren(:, deIDs) = [];
msh.elemGPx(:, deIDs)      = [];
msh.elemGPy(:, deIDs)      = [];
msh.elemGPz(:, deIDs)      = [];
msh.elemJac(:, deIDs)      = [];
msh.elemJxix(:, deIDs)     = [];
msh.elemJxiy(:, deIDs)     = [];
msh.elemJxiz(:, deIDs)     = [];
msh.elemJetax(:, deIDs)    = [];
msh.elemJetay(:, deIDs)    = [];
msh.elemJetaz(:, deIDs)    = [];
msh.elemJzetax(:, deIDs)   = [];
msh.elemJzetay(:, deIDs)   = [];
msh.elemJzetaz(:, deIDs)   = [];
msh.elemIsCurved(deIDs)    = [];
msh.faceNormalx(:, dfIDs)  = [];
msh.faceNormaly(:, dfIDs)  = [];
msh.faceNormalz(:, dfIDs)  = [];
msh.faceSize(dfIDs)        = [];
msh.faceType(dfIDs)        = [];
msh.faceNodes(:, dfIDs)    = [];
msh.faceEdges(:, dfIDs)    = [];
msh.faceElems(:, dfIDs)    = [];
msh.faceNums(:, dfIDs)     = [];
msh.faceR2L(dfIDs)         = [];
msh.faceLevel(dfIDs)       = []; 
msh.faceParent(dfIDs)      = []; 
msh.faceChildren(:, dfIDs) = []; 
msh.faceJac(:, dfIDs)      = [];
msh.faceIsCurved(dfIDs)    = [];
msh.edgeType(dedIDs)       = []; 
msh.edgeNodes(:, dedIDs)   = []; 
msh.edgeNElems(dedIDs)     = []; 
msh.edgeParent(dedIDs)     = []; 
msh.edgeChildren(:, dedIDs)= [];
msh.nodeCoor(:, dnIDs)     = [];
msh.massMatInv(deIDs)      = [];
msh.nElems = length(remEIDs);
msh.nFaces = length(remFIDs);
msh.nEdges = length(remEDIDs);
msh.nNodes = length(remNIDs);

%**************************************************************************
% Finally update the IDs of leaf elements and faces and get the new
% solution
%**************************************************************************
% Update the IDs of leaf elements
msh.LElems  = find(msh.elemChildren(1, :) == 0);
msh.nLElems = length(msh.LElems);

msh.elemLID             = zeros(1, msh.nElems);
msh.elemLID(msh.LElems) = 1 : msh.nLElems;

% Update the IDs of leaf faces
msh.LFaces  = find(msh.faceChildren(1, :) == 0);
msh.nLFaces = length(msh.LFaces);

msh.faceLID             = zeros(1, msh.nFaces);
msh.faceLID(msh.LFaces) = 1 : msh.nLFaces;

% Update the IDs of internal and various boundary faces
isInternal     = msh.faceType(msh.LFaces) == 0;
msh.intLFaces  = msh.LFaces(isInternal);
msh.nIntLFaces = length(msh.intLFaces);
bndLFaces      = msh.LFaces(~isInternal);
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = bndLFaces(msh.faceType(bndLFaces) == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end

% Update the solution in the new leaf elements
UT(:, (-nv + 1 : 0)' + nv * deIDs) = [];
U = UT(:, (-nv + 1 : 0)' + nv * msh.LElems);
if (layout == 1)
    U = U(:, (1 : nv : end)' + (0 : nv - 1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Inlined  Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine an edge
function [nodeID, nodeCoor, edgeChildren] = refineEdge(edgeID, edgeNodeCoors, addElement)

edgeChildren = msh.edgeChildren(:, edgeID);
if (edgeChildren(1) == 0)   
    edgeType     = msh.edgeType(edgeID);
    edgeNodes    = msh.edgeNodes(:, edgeID);
    edgeChildren = [msh.nEdges + 1; msh.nEdges + 2];
    
    nodeID = msh.nNodes + 1;
    nodeCoor = getEdgeNode(edgeNodeCoors(:, 1), edgeNodeCoors(:, 2), edgeType);
    msh.nodeCoor(:, nodeID) = nodeCoor;    
    msh.nNodes = msh.nNodes + 1;  
    
    msh.edgeChildren(:, edgeID)       = edgeChildren;    
    msh.edgeType(edgeChildren)        = edgeType;
    msh.edgeNodes(:, edgeChildren)    = [edgeNodes(1), nodeID; nodeID, edgeNodes(2)];
    msh.edgeNElems(edgeChildren)      = 1;
    msh.edgeParent(edgeChildren)      = edgeID;  
    msh.edgeChildren(:, edgeChildren) = 0;
    msh.nEdges = msh.nEdges + 2;
else
    nodeID   = msh.edgeNodes(2, edgeChildren(1));
    nodeCoor = msh.nodeCoor(:, nodeID);
    if addElement
        msh.edgeNElems(edgeChildren) = msh.edgeNElems(edgeChildren) + 1;
    end
end

end

% refine a face
function [elemNode, elemNodeCoor, elemEdges, faceIDs] = refineFace(elemID, faceID, elemEdges1, elemEdges2, elemNodes, subelemIDs, edgeTypes, elemNodeCoors4)
    
faceType     = msh.faceType(faceID);
faceElem     = msh.faceElems(:, faceID);
faceNum      = msh.faceNums(:, faceID);
faceR2L      = msh.faceR2L(faceID);
faceChildren = msh.faceChildren(:, faceID);
if (faceR2L == 1)
    if (faceChildren(1) == 0)
        faceChildren = msh.nFaces + 1 : msh.nFaces + 4;
        msh.faceChildren(:, faceID) = faceChildren';
        
        elemNode     = msh.nNodes + 1;
        elemNodes(5) = elemNode;
        msh.nNodes   = msh.nNodes + 1;
   
        elemEdges        = msh.nEdges + 1 : msh.nEdges + 4;
        elemEdges1(2, :) = elemEdges(1 : 2);
        elemEdges2(:, 2) = elemEdges(3 : 4)';
        if all(edgeTypes(1) == edgeTypes(2 : 4))
             msh.edgeType(elemEdges) = edgeTypes(1);
        else
            msh.edgeType(elemEdges) = 0;
        end        
        msh.edgeNodes(:, elemEdges)    = [elemNodes(2), elemNodes(5), elemNodes(4), elemNodes(5); ...
                                          elemNodes(5), elemNodes(8), elemNodes(5), elemNodes(6)];     
        msh.edgeNElems(elemEdges)      = 2;
        msh.edgeParent(elemEdges)      = 0;
        msh.edgeChildren(:, elemEdges) = 0;
        msh.nEdges = msh.nEdges + 4; 
        
        faceIDs = faceChildren;
        if (faceType == 0)  
            elemNodeCoor = 0.25 * sum(elemNodeCoors4, 2);
            faceNodes  = elemNodes;
            faceEdges1 = elemEdges1;
            faceEdges2 = elemEdges2;
        elseif (faceType == 1) 
            elemNodeCoor = 0.25 * sum(elemNodeCoors4, 2);        
            faceNodes4 = msh.faceNodes(:, faceID);
            if all(faceNodes4 == elemNodes([1, 3, 9, 7])')
                faceNodes  = elemNodes;
                faceEdges1 = elemEdges1;
                faceEdges2 = elemEdges2;
            else
                faceEdges4     = msh.faceEdges(:, faceID);
                faceNodeCoors4 = msh.nodeCoor(:, faceNodes4);
                faceEdgeNodes  = msh.edgeNodes(:, faceEdges4);
                faceNodes      = zeros(1, 9);
                faceEdges1     = zeros(3, 2);
                faceEdges2     = zeros(2, 3);
                
                faceNodes([1, 3, 9, 7]) = faceNodes4;     
                [faceNodes(4), ~, faceEdges1(1, :)] = refineEdge(faceEdges4(1), faceNodeCoors4(:, [1, 4]), true);
                [faceNodes(6), ~, faceEdges1(3, :)] = refineEdge(faceEdges4(2), faceNodeCoors4(:, 2 : 3), true);
                [faceNodes(2), ~, faceEdges2(:, 1)] = refineEdge(faceEdges4(3), faceNodeCoors4(:, 1 : 2), true);
                [faceNodes(8), ~, faceEdges2(:, 3)] = refineEdge(faceEdges4(4), faceNodeCoors4(:, 3 : 4), true);
                if (faceEdgeNodes(1, 1) == faceNodes(7))
                    faceEdges1(1, :) = flip(faceEdges1(1, :));
                end
                if (faceEdgeNodes(1, 2) == faceNodes(9))
                    faceEdges1(3, :) = flip(faceEdges1(3, :));
                end   
                if (faceEdgeNodes(1, 3) == faceNodes(3))
                    faceEdges2(:, 1) = flip(faceEdges2(:, 1));
                end  
                if (faceEdgeNodes(1, 4) == faceNodes(9))
                    faceEdges2(:, 3) = flip(faceEdges2(:, 3));
                end 
                
                faceNodes(5) = msh.nNodes + 1;
                msh.nodeCoor(:, faceNodes(5)) = 0.25 * sum(msh.nodeCoor(:, faceNodes4), 2); 
                msh.nNodes = msh.nNodes + 1;
                
                faceEdges        = msh.nEdges + 1 : msh.nEdges + 4;
                faceEdges1(2, :) = faceEdges(1 : 2);
                faceEdges2(:, 2) = faceEdges(3 : 4)';
                msh.edgeType(faceEdges)        = 0;
                msh.edgeNodes(:, faceEdges)    = [faceNodes(2), faceNodes(5), faceNodes(4), faceNodes(5); ...
                                                  faceNodes(5), faceNodes(8), faceNodes(5), faceNodes(6)];     
                msh.edgeNElems(faceEdges)      = 2;
                msh.edgeParent(faceEdges)      = 0;
                msh.edgeChildren(:, faceEdges) = 0;
                msh.nEdges = msh.nEdges + 4; 
            end           
        else
            elemNodeCoor = getFaceNode(elemNodeCoors4(:, 1), elemNodeCoors4(:, 2), elemNodeCoors4(:, 3), elemNodeCoors4(:, 4), faceType);
            faceNodes  = elemNodes;
            faceEdges1 = elemEdges1;
            faceEdges2 = elemEdges2;           
        end
        msh.nodeCoor(:, elemNode) = elemNodeCoor;
        faceNodeCoors = msh.nodeCoor(:, faceNodes);
        
        msh.faceNodes(1, faceChildren) = faceNodes([1, 2, 4, 5]);
        msh.faceNodes(2, faceChildren) = faceNodes([2, 3, 5, 6]);
        msh.faceNodes(3, faceChildren) = faceNodes([5, 6, 8, 9]);
        msh.faceNodes(4, faceChildren) = faceNodes([4, 5, 7, 8]);
        msh.faceEdges(1, faceChildren) = reshape(faceEdges1(1 : 2, :), [1, 4]);
        msh.faceEdges(2, faceChildren) = reshape(faceEdges1(2 : 3, :), [1, 4]);
        msh.faceEdges(3, faceChildren) = reshape(faceEdges2(:, 1 : 2), [1, 4]);
        msh.faceEdges(4, faceChildren) = reshape(faceEdges2(:, 2 : 3), [1, 4]);
        if (faceType == 0 || faceType == 1)
            if (faceElem(1) == elemID)
                msh.faceElems(1, faceChildren) = subelemIDs;
                msh.faceElems(2, faceChildren) = faceElem(2);
                msh.faceNums( 1, faceChildren) = faceNum(1);
                msh.faceNums( 2, faceChildren) = 4 * faceNum(2) + 3 : 4 * faceNum(2) + 6;
            else
                msh.faceElems(1, faceChildren) = faceElem(1);
                msh.faceElems(2, faceChildren) = subelemIDs;
                msh.faceNums( 1, faceChildren) = 4 * faceNum(1) + 3 : 4 * faceNum(1) + 6;
                msh.faceNums( 2, faceChildren) = faceNum(2);
            end             
        else
            msh.faceElems(1, faceChildren) = subelemIDs;
            msh.faceNums( 1, faceChildren) = faceNum(1);
        end
        msh.faceR2L(faceChildren)  = faceR2L;
        
        tempNormals = 0.5 * mycross3(faceNodeCoors(:, [1, 2, 4, 5]) - faceNodeCoors(:, [5, 6, 8, 9]), faceNodeCoors(:,[2, 3, 5, 6]) - faceNodeCoors(:, [4, 5, 7, 8]));
        tempSizes   = vecnorm(tempNormals);
        tempNormals = tempNormals ./ tempSizes;
        if mod(faceNum(1), 2) == 1
            tempNormals = -tempNormals;
        end       
        msh.faceNormalx(faceChildren) = tempNormals(1, :);
        msh.faceNormaly(faceChildren) = tempNormals(2, :);
        msh.faceNormalz(faceChildren) = tempNormals(3, :);
        msh.faceSize(faceChildren)    = tempSizes;
        msh.faceType(faceChildren)    = faceType;
        msh.faceLevel(faceChildren)   = msh.faceLevel(faceID) + 1;
        msh.faceParent(faceChildren)  = faceID;
        msh.faceJac(faceChildren)     = tempSizes / 4;
        if msh.faceIsCurved(faceID)
            msh.faceIsCurved(faceChildren) = true;
            msh.elemIsCurved(subelemIDs) = true;
        end
        msh.nFaces = msh.nFaces + 4;
    else  
        % It must be an internal or periodic face in this case
        if (faceType == 0)  
            elemNode = msh.faceNodes(3, faceChildren(1));
            elemNodeCoor = msh.nodeCoor(:, elemNode);
            elemEdges = [msh.faceEdges(2, faceChildren(1)), msh.faceEdges(1, faceChildren(4)), msh.faceEdges(4, faceChildren(1)), msh.faceEdges(3, faceChildren(4))];
            msh.edgeNElems(elemEdges) = 4;
        elseif (faceType == 1)
            if all(msh.faceNodes(:, faceID) == elemNodes([1, 3, 9, 7])')
                elemNode = msh.faceNodes(3, faceChildren(1));
                elemNodeCoor = msh.nodeCoor(:, elemNode);
                elemEdges = [msh.faceEdges(2, faceChildren(1)), msh.faceEdges(1, faceChildren(4)), msh.faceEdges(4, faceChildren(1)), msh.faceEdges(3, faceChildren(4))];
                msh.edgeNElems(elemEdges) = 4;
            else
                edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
                msh.edgeNElems(edgeIDs8) = msh.edgeNElems(edgeIDs8) + 1;  

                elemNode = msh.nNodes + 1;
                elemNodeCoor = 0.25 * sum(elemNodeCoors4, 2);
                msh.nodeCoor(:, elemNode) = elemNodeCoor; 
                msh.nNodes = msh.nNodes + 1;

                elemEdges = msh.nEdges + 1 : msh.nEdges + 4;
                msh.edgeType(elemEdges)        = 0;
                msh.edgeNodes(:, elemEdges)    = [elemNodes(2),     elemNode, elemNodes(4),     elemNode; ...
                                                      elemNode, elemNodes(8),     elemNode, elemNodes(6)];     
                msh.edgeNElems(elemEdges)      = 2;
                msh.edgeParent(elemEdges)      = 0;
                msh.edgeChildren(:, elemEdges) = 0;
                msh.nEdges = msh.nEdges + 4;
            end
        else
            error('Wrong face type')
        end    
        faceIDs = faceChildren;
        
        if (faceElem(1) == elemID)
            msh.faceElems(1, faceChildren) = subelemIDs;
            msh.faceNums( 1, faceChildren) = faceNum(1);
        else
            msh.faceElems(2, faceChildren) = subelemIDs;
            msh.faceNums( 2, faceChildren) = faceNum(2);               
        end        
    end
else
    faceNodes4 = msh.faceNodes(:, faceID);
    if (faceChildren(1) == 0)
        faceChildren = msh.nFaces + 1 : msh.nFaces + 4;
        msh.faceChildren(:, faceID) = faceChildren';
        
        elemNode     = msh.nNodes + 1;
        elemNodes(5) = elemNode;
        msh.nNodes   = msh.nNodes + 1;   
        
        elemEdges        = msh.nEdges + 1 : msh.nEdges + 4;
        elemEdges1(2, :) = elemEdges(1 : 2);
        elemEdges2(:, 2) = elemEdges(3 : 4)';
        if all(edgeTypes(1) == edgeTypes(2 : 4))
             msh.edgeType(elemEdges) = edgeTypes(1);
        else
            msh.edgeType(elemEdges) = 0;
        end        
        msh.edgeNodes(:, elemEdges)    = [elemNodes(2), elemNodes(5), elemNodes(4), elemNodes(5); ...
                                          elemNodes(5), elemNodes(8), elemNodes(5), elemNodes(6)];     
        msh.edgeNElems(elemEdges)      = 2;
        msh.edgeParent(elemEdges)      = 0;
        msh.edgeChildren(:, elemEdges) = 0;        
        msh.nEdges = msh.nEdges + 4;  
        
        if (faceElem(1) == elemID) && all(faceNodes4 == elemNodes([1, 3, 9, 7])')
            faceNodes  = elemNodes;
            faceEdges1 = elemEdges1;
            faceEdges2 = elemEdges2;        
        elseif (faceElem(2) == elemID) && all(faceNodes4 == reorderNodeIDs1_r2l(elemNodes([1, 3, 9, 7])', faceR2L))
            faceNodes  = reorderNodeIDs2_r2l(elemNodes, faceR2L);
            [faceEdges1, faceEdges2] = reorderEdgeIDs1_r2l(elemEdges1, elemEdges2, faceR2L);
        else
            faceEdges4     = msh.faceEdges(:, faceID);
            faceNodeCoors4 = msh.nodeCoor(:, faceNodes4);
            faceEdgeNodes  = msh.edgeNodes(:, faceEdges4);
            faceNodes      = zeros(1, 9);
            faceEdges1     = zeros(3, 2);
            faceEdges2     = zeros(2, 3);

            faceNodes([1, 3, 9, 7]) = faceNodes4;     
            [faceNodes(4), ~, faceEdges1(1, :)] = refineEdge(faceEdges4(1), faceNodeCoors4(:, [1, 4]), true);
            [faceNodes(6), ~, faceEdges1(3, :)] = refineEdge(faceEdges4(2), faceNodeCoors4(:, 2 : 3), true);
            [faceNodes(2), ~, faceEdges2(:, 1)] = refineEdge(faceEdges4(3), faceNodeCoors4(:, 1 : 2), true);
            [faceNodes(8), ~, faceEdges2(:, 3)] = refineEdge(faceEdges4(4), faceNodeCoors4(:, 3 : 4), true);
            if (faceEdgeNodes(1, 1) == faceNodes(7))
                faceEdges1(1, :) = flip(faceEdges1(1, :));
            end
            if (faceEdgeNodes(1, 2) == faceNodes(9))
                faceEdges1(3, :) = flip(faceEdges1(3, :));
            end   
            if (faceEdgeNodes(1, 3) == faceNodes(3))
                faceEdges2(:, 1) = flip(faceEdges2(:, 1));
            end  
            if (faceEdgeNodes(1, 4) == faceNodes(9))
                faceEdges2(:, 3) = flip(faceEdges2(:, 3));
            end 

            faceNodes(5) = msh.nNodes + 1;
            msh.nodeCoor(:, faceNodes(5)) = 0.25 * sum(msh.nodeCoor(:, faceNodes4), 2); 
            msh.nNodes = msh.nNodes + 1;

            faceEdges        = msh.nEdges + 1 : msh.nEdges + 4;
            faceEdges1(2, :) = faceEdges(1 : 2);
            faceEdges2(:, 2) = faceEdges(3 : 4)';
            msh.edgeType(faceEdges)        = 0;
            msh.edgeNodes(:, faceEdges)    = [faceNodes(2), faceNodes(5), faceNodes(4), faceNodes(5); ...
                                              faceNodes(5), faceNodes(8), faceNodes(5), faceNodes(6)];     
            msh.edgeNElems(faceEdges)      = 2;
            msh.edgeParent(faceEdges)      = 0;
            msh.edgeChildren(:, faceEdges) = 0;
            msh.nEdges = msh.nEdges + 4;                 
        end         
        
        msh.faceNodes(1, faceChildren) = faceNodes([1, 2, 4, 5]);
        msh.faceNodes(2, faceChildren) = faceNodes([2, 3, 5, 6]);
        msh.faceNodes(3, faceChildren) = faceNodes([5, 6, 8, 9]);
        msh.faceNodes(4, faceChildren) = faceNodes([4, 5, 7, 8]);
        msh.faceEdges(1, faceChildren) = reshape(faceEdges1(1 : 2, :), [1, 4]);
        msh.faceEdges(2, faceChildren) = reshape(faceEdges1(2 : 3, :), [1, 4]);
        msh.faceEdges(3, faceChildren) = reshape(faceEdges2(:, 1 : 2), [1, 4]);
        msh.faceEdges(4, faceChildren) = reshape(faceEdges2(:, 2 : 3), [1, 4]);        
        if (faceType == 0 || faceType == 1)
            elemNodeCoor = 0.25 * sum(elemNodeCoors4, 2);            
            if (faceElem(1) == elemID)
                faceIDs = faceChildren;
                msh.faceElems(1, faceChildren) = subelemIDs;
                msh.faceElems(2, faceChildren) = faceElem(2);
                msh.faceNums( 1, faceChildren) = faceNum(1);
                msh.faceNums( 2, faceChildren) = reorderFaceIDs_r2l(4 * faceNum(2) + 3 : 4 * faceNum(2) + 6, faceR2L);                
            else
                faceIDs = reorderFaceIDs_l2r(faceChildren, faceR2L);
                msh.faceElems(1, faceChildren) = faceElem(1);
                msh.faceElems(2, faceChildren) = reorderFaceIDs_r2l(subelemIDs, faceR2L);
                msh.faceNums( 1, faceChildren) = 4 * faceNum(1) + 3 : 4 * faceNum(1) + 6;
                msh.faceNums( 2, faceChildren) = faceNum(2);                
            end
        else
            elemNodeCoor = getFaceNode(elemNodeCoors4(:, 1), elemNodeCoors4(:, 2), elemNodeCoors4(:, 3), elemNodeCoors4(:, 4), faceType);
            faceIDs = faceChildren;
            msh.faceElems(1, faceChildren) = subelemIDs;
            msh.faceNums( 1, faceChildren) = faceNum(1);        
        end      
        msh.faceR2L(faceChildren)  = faceR2L;

        msh.nodeCoor(:, elemNode) = elemNodeCoor; 
        faceNodeCoors = msh.nodeCoor(:, faceNodes);
        tempNormals = 0.5 * mycross3(faceNodeCoors(:, [1, 2, 4, 5]) - faceNodeCoors(:, [5, 6, 8, 9]), faceNodeCoors(:,[2, 3, 5, 6]) - faceNodeCoors(:, [4, 5, 7, 8]));
        tempSizes   = vecnorm(tempNormals);
        tempNormals = tempNormals ./ tempSizes;
        if mod(faceNum(1), 2) == 1
            tempNormals = -tempNormals;
        end       
        msh.faceNormalx(faceChildren) = tempNormals(1, :);
        msh.faceNormaly(faceChildren) = tempNormals(2, :);
        msh.faceNormalz(faceChildren) = tempNormals(3, :);
        msh.faceSize(faceChildren)    = tempSizes;
        msh.faceType(faceChildren)    = faceType;
        msh.faceLevel(faceChildren)   = msh.faceLevel(faceID) + 1;
        msh.faceParent(faceChildren)  = faceID;
        msh.faceJac(faceChildren)     = tempSizes / 4;
        if msh.faceIsCurved(faceID)
            msh.faceIsCurved(faceChildren) = true;
            msh.elemIsCurved(subelemIDs) = true;
        end
        msh.nFaces = msh.nFaces + 4;        
    else
        % It must be an internal or periodic face in this case
        if (faceElem(1) == elemID) && all(faceNodes4 == elemNodes([1, 3, 9, 7])')
            elemNode = msh.faceNodes(3, faceChildren(1));
            elemNodeCoor = msh.nodeCoor(:, elemNode);
            elemEdges = [msh.faceEdges(2, faceChildren(1)), msh.faceEdges(1, faceChildren(4)), msh.faceEdges(4, faceChildren(1)), msh.faceEdges(3, faceChildren(4))];
            msh.edgeNElems(elemEdges) = 4;
        elseif (faceElem(2) == elemID) && all(faceNodes4 == reorderNodeIDs1_r2l(elemNodes([1, 3, 9, 7])', faceR2L))
            elemNode = msh.faceNodes(3, faceChildren(1));
            elemNodeCoor = msh.nodeCoor(:, elemNode);
            elemEdges = reorderEdgeIDs2_l2r([msh.faceEdges(2, faceChildren(1)), msh.faceEdges(1, faceChildren(4)), msh.faceEdges(4, faceChildren(1)), msh.faceEdges(3, faceChildren(4))], faceR2L);
            msh.edgeNElems(elemEdges) = 4;
        else
            edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
            msh.edgeNElems(edgeIDs8) = msh.edgeNElems(edgeIDs8) + 1;  
            
            elemNode = msh.nNodes + 1;
            elemNodeCoor = 0.25 * sum(elemNodeCoors4, 2); 
            msh.nodeCoor(:, elemNode) = elemNodeCoor; 
            msh.nNodes = msh.nNodes + 1;
            
            elemEdges = msh.nEdges + 1 : msh.nEdges + 4;
            msh.edgeType(elemEdges)        = 0;
            msh.edgeNodes(:, elemEdges)    = [elemNodes(2),     elemNode, elemNodes(4),     elemNode; ...
                                                  elemNode, elemNodes(8),     elemNode, elemNodes(6)];     
            msh.edgeNElems(elemEdges)      = 2;
            msh.edgeParent(elemEdges)      = 0;
            msh.edgeChildren(:, elemEdges) = 0;
            msh.nEdges = msh.nEdges + 4;
        end
        
        if (faceElem(1) == elemID)
            faceIDs = faceChildren;
            msh.faceElems(1, faceChildren) = subelemIDs;
            msh.faceNums( 1, faceChildren) = faceNum(1);
        else
            faceIDs = reorderFaceIDs_l2r(faceChildren, faceR2L);
            msh.faceElems(2, faceChildren) = reorderFaceIDs_r2l(subelemIDs, faceR2L);
            msh.faceNums( 2, faceChildren) = faceNum(2);                
        end                
    end
end

end

% coarsen a face
function coarsenFace(elemID, faceID, edgeIDs, nodeIDs)
    
faceType     = msh.faceType(faceID);
faceElem     = msh.faceElems(:, faceID);
faceNum      = msh.faceNums(:, faceID);
faceR2L      = msh.faceR2L(faceID);
faceChildren = msh.faceChildren(:, faceID);   
if (faceType == 0 || faceType == 1)
    isLeft = faceElem(1) == elemID;
    adjElemID = faceElem(1) * ~isLeft + faceElem(2) * isLeft;   
    if (msh.elemChildren(1, adjElemID) ~= 0)
        if isLeft
            msh.faceElems(1, faceChildren) = elemID;
            msh.faceNums( 1, faceChildren) = 4 * faceNum(1) + 3 : 4 * faceNum(1) + 6;
            if (faceType == 0)
                msh.edgeNElems(edgeIDs) = 2;
            else
                if all(msh.faceNodes(:, faceID) == nodeIDs(1 : 4)')
                    msh.edgeNElems(edgeIDs) = 2;
                else
                    edgeIDs4 = [msh.faceEdges([2, 4], faceChildren(1)); msh.faceEdges([1, 3], faceChildren(4))];
                    msh.edgeNElems(edgeIDs4) = 2;
                    edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
                    msh.edgeNElems(edgeIDs8) = msh.edgeNElems(edgeIDs8) - 1;

                    dedIDs(nded + 1 : nded + 4) = edgeIDs;
                    nded = nded + 4;

                    dnIDs(ndn + 1) = nodeIDs(5);
                    ndn = ndn + 1;               
                end 
            end
        else
            msh.faceElems(2, faceChildren) = elemID;
            msh.faceNums( 2, faceChildren) = reorderFaceIDs_r2l(4 * faceNum(2) + 3 : 4 * faceNum(2) + 6, faceR2L);
            if (faceType == 0)
                msh.edgeNElems(edgeIDs) = 2;
            else
                if all(msh.faceNodes(:, faceID) == reorderNodeIDs1_r2l(nodeIDs(1 : 4)', faceR2L))
                    msh.edgeNElems(edgeIDs) = 2;
                else
                    edgeIDs4 = [msh.faceEdges([2, 4], faceChildren(1)); msh.faceEdges([1, 3], faceChildren(4))];
                    msh.edgeNElems(edgeIDs4) = 2;
                    edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
                    msh.edgeNElems(edgeIDs8) = msh.edges(3, edgeIDs8) - 1; 

                    dedIDs(nded + 1 : nded + 4) = edgeIDs;
                    nded = nded + 4;

                    dnIDs(ndn + 1) = nodeIDs(5);
                    ndn = ndn + 1;                       
                end
            end
        end
    else
        msh.faceChildren(:, faceID) = 0;
        dfIDs(ndf + 1 : ndf + 4) = faceChildren;
        ndf = ndf + 4;
        
        dedIDs(nded + 1 : nded + 4) = edgeIDs;
        nded = nded + 4;
        
        dnIDs(ndn + 1) = nodeIDs(5);
        ndn = ndn + 1;
        if (faceType == 1)
            if isLeft
                if ~all(msh.faceNodes(:, faceID) == nodeIDs(1 : 4)')
                    dedIDs(nded + 1 : nded + 4) = [msh.faceEdges([2, 4], faceChildren(1)); msh.faceEdges([1, 3], faceChildren(4))];
                    nded = nded + 4;

                    dnIDs(ndn + 1) = msh.faceNodes(3, faceChildren(1));
                    ndn = ndn + 1;                

                    edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
                    nae_ = msh.edgeNElems(edgeIDs8) - 1;
                    msh.edgeNElems(edgeIDs8) = nae_;  

                    delEdgeIDs8 = edgeIDs8(:, nae_(1 : 2 : end) == 0);
                    if ~isempty(delEdgeIDs8)
                        nded8 = numel(delEdgeIDs8);
                        dedIDs(nded + 1 : nded + nded8) = delEdgeIDs8;
                        nded = nded + nded8;

                        dnIDs(ndn + 1 : ndn + 0.5 * nded8) = msh.edgeNodes(2, delEdgeIDs8(1, :));
                        ndn = ndn + 0.5 * nded8;
                    end  
                end
            else
                if ~all(msh.faceNodes(:, faceID) == reorderNodeIDs1_r2l(nodeIDs(1 : 4)', faceR2L))
                    dedIDs(nded + 1 : nded + 4) = [msh.faceEdges([2, 4], faceChildren(1)); msh.faceEdges([1, 3], faceChildren(4))];
                    nded = nded + 4;

                    dnIDs(ndn + 1) = msh.faceNodes(3, faceChildren(1));
                    ndn = ndn + 1;

                    edgeIDs8 = msh.edgeChildren(:, msh.faceEdges(:, faceID));
                    nae_ = msh.edgeNElems(edgeIDs8) - 1;
                    msh.edgeNElems(edgeIDs8) = nae_;  

                    delEdgeIDs8 = edgeIDs8(:, nae_(1 : 2 : end) == 0);
                    if ~isempty(delEdgeIDs8)
                        nded8 = numel(delEdgeIDs8);
                        dedIDs(nded + 1 : nded + nded8) = delEdgeIDs8;
                        nded = nded + nded8;

                        dnIDs(ndn + 1 : ndn + 0.5 * nded8) = msh.edgeNodes(2, delEdgeIDs8(1, :));
                        ndn = ndn + 0.5 * nded8;
                    end                   
                end          
            end
        end
    end
else
    msh.faceChildren(:, faceID) = 0;
    dfIDs(ndf + 1 : ndf + 4) = faceChildren;
    ndf = ndf + 4;

    dedIDs(nded + 1 : nded + 4) = edgeIDs;
    nded = nded + 4;

    dnIDs(ndn + 1) = nodeIDs(5);
    ndn = ndn + 1;    
end    

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given the subfaces IDs of one face which associate with the left element, 
% and the face-to-face index, reorder them according to the right element
function IDs = reorderFaceIDs_l2r(IDs, f2f)

switch f2f
    case 1
    case 2
        IDs = IDs([2, 4, 1, 3]);
    case 3
        IDs = IDs([4, 3, 2, 1]);
    case 4
        IDs = IDs([3, 1, 4, 2]);
    case 5
        IDs = IDs([2, 1, 4, 3]);
    case 6
        IDs = IDs([3, 4, 1, 2]);
    case 7
        IDs = IDs([1, 3, 2, 4]);
    case 8
        IDs = IDs([4, 2, 3, 1]);
    otherwise
        error('Wrong face to face index')
end    

end

% Given the subfaces IDs of one face which associate with the right element, 
% and the face-to-face index, reorder them according to the left element
function IDs = reorderFaceIDs_r2l(IDs, f2f)
    
switch f2f
    case 1
    case 2
        IDs = IDs([3, 1, 4, 2]);
    case 3
        IDs = IDs([4, 3, 2, 1]);
    case 4
        IDs = IDs([2, 4, 1, 3]);
    case 5
        IDs = IDs([2, 1, 4, 3]);
    case 6
        IDs = IDs([3, 4, 1, 2]);
    case 7
        IDs = IDs([1, 3, 2, 4]);
    case 8
        IDs = IDs([4, 2, 3, 1]);
    otherwise
        error('Wrong face to face index')
end

end

% Given the four face node IDs ordered w.r.t. the rifght element, and the
% face-to-face index, reorder them according to the left element
function IDs = reorderNodeIDs1_r2l(IDs, f2f)
    
switch f2f
    case 1
    case 2
        IDs = IDs([4, 1, 2, 3]);
    case 3
        IDs = IDs([3, 4, 1, 2]);
    case 4
        IDs = IDs([2, 3, 4, 1]);
    case 5
        IDs = IDs([2, 1, 4, 3]);
    case 6
        IDs = IDs([4, 3, 2, 1]);
    case 7
        IDs = IDs([1, 4, 3, 2]);
    case 8
        IDs = IDs([3, 2, 1, 4]);
    otherwise
        error('Wrong face to face index')
end

end

% Given the nine face node IDs ordered w.r.t. the rifght element, and the
% face-to-face index, reorder them according to the left element
function IDs = reorderNodeIDs2_r2l(IDs, f2f)
    
si = size(IDs);
IDs = reshape(IDs, [3, 3]);    
switch f2f
    case 1
    case 2
        IDs = rot90(IDs, 1);
    case 3
        IDs = rot90(IDs, 2);
    case 4
        IDs = rot90(IDs, 3);
    case 5
        IDs = flipud(IDs);
    case 6
        IDs = fliplr(IDs);
    case 7
        IDs = IDs';
    case 8
        IDs = rot90(IDs', 2);
    otherwise
        error('Wrong face to face index')
end
IDs = reshape(IDs, si);

end

% Given the face edge IDs ordered w.r.t. the rifght element, and the
% face-to-face index, reorder them according to the left element
function [IDs1, IDs2] = reorderEdgeIDs1_r2l(IDs1_, IDs2_, f2f)
    
switch f2f
    case 1
    case 2
        IDs1 = rot90(IDs2_);
        IDs2 = rot90(IDs1_);
    case 3
        IDs1 = rot90(IDs1_, 2);
        IDs2 = rot90(IDs2_, 2);
    case 4
        IDs1 = rot90(IDs2_, 3);
        IDs2 = rot90(IDs1_, 3);
    case 5
        IDs1 = flipud(IDs1_);
        IDs2 = flipud(IDs2_);
    case 6
        IDs1 = fliplr(IDs1_);
        IDs2 = fliplr(IDs2_);
    case 7
        IDs1 = IDs2_';
        IDs2 = IDs1_';
    case 8
        IDs1 = rot90(IDs2_', 2);
        IDs2 = rot90(IDs1_', 2);
    otherwise
        error('Wrong face to face index')
end

end

% Given the internal face edge IDs ordered w.r.t. the left element, and the
% face-to-face index, reorder them according to the right element
function IDs = reorderEdgeIDs2_l2r(IDs, f2f)
    
switch f2f
    case 1
    case 2
        IDs = IDs([4, 3, 1, 2]);
    case 3
        IDs = IDs([2, 1, 4, 3]);
    case 4
        IDs = IDs([3, 4, 2, 1]);
    case 5
        IDs = IDs([1, 2, 4, 3]);
    case 6
        IDs = IDs([2, 1, 3, 4]);
    case 7
        IDs = IDs([3, 4, 1, 2]);
    case 8
        IDs = IDs([4, 3, 2, 1]);
    otherwise
        error('Wrong face to face index')
end

end












