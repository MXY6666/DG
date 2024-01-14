% Suppose the level difference between neighboring elements is no more than
% 1 when refining and coarsening
% Copyright (C) Xiangyi Meng.

function [msh, U] = refineAndCoarsenQuadMesh(msh, candReIDs, candCeIDs, nv, U, quad, bs, getBndNode, layout)

if (nargin < 1)
    error('Not enough arguments')
end

if (msh.type ~= 212)
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
    U = zeros(3, msh.nLElems);
end
if (size(U, 2) ~= msh.nLElems * nv)
    error('The given mesh and the solution coefficients does not coinside')
end

if (nargin < 6) || isempty(quad)
    quad = GaussQuadratureRule_quad([2, 2], 201);
end

if (nargin < 7) || isempty(bs)
    quad1 = GaussQuadratureRule_line(2, 101);
    bs = setBasisFunctionSet_square(quad1, quad, 201, 2);
    bs = addBasisFunctionSet_square(bs, quad1, quad);
end
if (size(U, 1) ~= bs.nb)
    error('The given U and bs does not coinside')
end

if (nargin < 8) || isempty(getBndNode)
    getBndNode = @(p1, p2, faceType)0.5 * (p1 + p2);
end

if (nargin < 9) || isempty(layout)
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
    elemFaces = msh.elemFaces(:, ID);   
    faceTypes = msh.faceType(elemFaces);   
    faceElems = msh.faceElems(:, elemFaces);
    
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

% Then we refine the elements by adding new elements, faces and nodes
% preallocation for elements
msh.elemSize(msh.nElems + 1 : msh.nElems + 4 * nre)        = 0;
msh.elemDiam(msh.nElems + 1 : msh.nElems + 4 * nre)        = 0;
msh.elemNodes(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0;
msh.elemFaces(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0;
msh.elemLevel(msh.nElems + 1 : msh.nElems + 4 * nre)       = 0;
msh.elemParent(msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemChildren(:, msh.nElems + 1 : msh.nElems + 4 * nre) = 0;
msh.elemGPx(:, msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemGPy(:, msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemJac(:, msh.nElems + 1 : msh.nElems + 4 * nre)      = 0;
msh.elemJxix(:, msh.nElems + 1 : msh.nElems + 4 * nre)     = 0; 
msh.elemJxiy(:, msh.nElems + 1 : msh.nElems + 4 * nre)     = 0; 
msh.elemJetax(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0; 
msh.elemJetay(:, msh.nElems + 1 : msh.nElems + 4 * nre)    = 0; 
msh.elemIsCurved(msh.nElems + 1 : msh.nElems + 4 * nre)    = false;
% preallocation for faces
msh.faceNormalx(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)  = 0; 
msh.faceNormaly(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)  = 0; 
msh.faceSize(msh.nFaces + 1 : msh.nFaces + 12 * nre)        = 0; 
msh.faceType(msh.nFaces + 1 : msh.nFaces + 12 * nre)        = 0;
msh.faceNodes(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)    = 0;
msh.faceElems(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)    = 0; 
msh.faceNums(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)     = 0; 
msh.faceLevel(msh.nFaces + 1 : msh.nFaces + 12 * nre)       = 0; 
msh.faceParent(msh.nFaces + 1 : msh.nFaces + 12 * nre)      = 0; 
msh.faceChildren(:, msh.nFaces + 1 : msh.nFaces + 12 * nre) = 0; 
msh.faceJac(:, msh.nFaces + 1 : msh.nFaces + 12 * nre)      = 0; 
msh.faceIsCurved(msh.nFaces + 1 : msh.nFaces + 12 * nre)    = false; 
% preallocation for ndoes
msh.nodeCoor(:, msh.nNodes + 1 : msh.nNodes + 6 * nre) = 0; 
% preallocation for inverse of mass matrices
msh.massMatInv{msh.nElems + 1 + 4 * nre} = [];
% preallocation for solution
UT = zeros(size(U, 1), (msh.nElems + 4 * nre) * nv);
if (layout == 1)
    U = U(:, (0 : nv - 1)' * msh.nLElems + (1 : msh.nLElems));
end
UT(:, (-nv + 1 : 0)' + nv * msh.LElems) = U;

subfIDs   = zeros(2, 6);
nodeIDs   = zeros(1, 5);
nodeCoors = zeros(2, 5);
for i = 1 : nre
    % current loop element and its nodes and faces
    ID = reIDs(i);
    elemNodes     = msh.elemNodes(:, ID);
    elemFaces     = msh.elemFaces(:, ID);    
    elemNodeCoors = msh.nodeCoor(:, elemNodes);
    faceNormalxs  = msh.faceNormalx(elemFaces);
    faceNormalys  = msh.faceNormaly(elemFaces);
    faceSizes     = msh.faceSize(elemFaces);
    faceTypes     = msh.faceType(elemFaces);
    faceNodes     = msh.faceNodes(:, elemFaces);
    faceElems     = msh.faceElems(:, elemFaces);
    faceNums      = msh.faceNums(:, elemFaces);
    faceLevels    = msh.faceLevel(elemFaces);
    faceChildren  = msh.faceChildren(:, elemFaces);  
    faceIsCurved  = msh.faceIsCurved(elemFaces);
       
    % creation of subelements
    elemChildren = msh.nElems + 1 : msh.nElems + 4;
    msh.elemChildren(:, ID)      = elemChildren';
    msh.elemLevel(elemChildren)  = msh.elemLevel(ID) + 1;
    msh.elemParent(elemChildren) = ID;

    % creation of subfaces and new nodes
    % first face
    if (faceChildren(1, 1) == 0)
        faceChildren(:, 1) = [msh.nFaces + 1; msh.nFaces + 2];
        msh.faceChildren(:, elemFaces(1)) = faceChildren(:, 1);
        
        nodeIDs(1) = msh.nNodes + 1;
        msh.nNodes = msh.nNodes + 1;
        if (faceTypes(1) == 0)
            nodeCoors(:, 1) = 0.5 * sum(elemNodeCoors(:, [4, 1]), 2);
            nodeID = nodeIDs(1);
        elseif (faceTypes(1) == 1)
            nodeCoors(:, 1) = 0.5 * sum(elemNodeCoors(:, [4, 1]), 2);
            if (faceNodes(1, 1) == elemNodes(4) && faceNodes(2, 1) == elemNodes(1))
                nodeID = nodeIDs(1);
            elseif (faceNodes(1, 1) == elemNodes(1) && faceNodes(2, 1) == elemNodes(4))
                nodeID = nodeIDs(1);
            else
                nodeID = msh.nNodes + 1;
                msh.nodeCoor(:, nodeID) = 0.5 * sum(msh.nodeCoor(:, faceNodes(:, 1)), 2);
                msh.nNodes = msh.nNodes + 1; 
            end
        else
            nodeCoors(:, 1) = getBndNode(elemNodeCoors(:, 4), elemNodeCoors(:, 1), faceTypes(1));
            nodeID = nodeIDs(1);
        end
        msh.nodeCoor(:, nodeIDs(1)) = nodeCoors(:, 1);
        
        msh.faceNodes(:, faceChildren(:, 1)) = [faceNodes(1, 1), nodeID; nodeID, faceNodes(2, 1)];        
        if (faceTypes(1) == 0 || faceTypes(1) == 1)        
            if (faceElems(1, 1) == ID)          
                subfIDs(:, 1) = faceChildren(:, 1);
                msh.faceElems(1, faceChildren(:, 1)) = elemChildren([3, 1]);
                msh.faceElems(2, faceChildren(:, 1)) = faceElems(2, 1);
                msh.faceNums( 1, faceChildren(:, 1)) = 1;
                msh.faceNums( 2, faceChildren(:, 1)) = [2 * faceNums(2, 1) + 4, 2 * faceNums(2, 1) + 3]; 
            else              
                subfIDs(:, 1) = faceChildren([2, 1], 1);
                msh.faceElems(1, faceChildren(:, 1)) = faceElems(1, 1);
                msh.faceElems(2, faceChildren(:, 1)) = elemChildren([1, 3]);
                msh.faceNums( 1, faceChildren(:, 1)) = [2 * faceNums(1, 1) + 3, 2 * faceNums(1, 1) + 4];
                msh.faceNums( 2, faceChildren(:, 1)) = 1;
            end           
        else  
            subfIDs(:, 1) = faceChildren(:, 1);
            msh.faceElems(1, faceChildren(:, 1)) = elemChildren([3, 1]);
            msh.faceNums( 1, faceChildren(:, 1)) = 1;           
        end     
        if faceIsCurved(1)
            tempNormals = [nodeCoors(2, 1) - elemNodeCoors(2, 4), elemNodeCoors(2, 1) - nodeCoors(2, 1); elemNodeCoors(1, 4) - nodeCoors(1, 1), nodeCoors(1, 1) - elemNodeCoors(1, 1)];
            tempSizes   = vecnorm(tempNormals);
            tempNormals = tempNormals ./ tempSizes;
            
            msh.faceNormalx(faceChildren(:, 1))  = tempNormals(1, :);
            msh.faceNormaly(faceChildren(:, 1))  = tempNormals(2, :);
            msh.faceSize(faceChildren(:, 1))     = tempSizes;
            msh.faceJac(faceChildren(:, 1))      = tempSizes / 2; 
            msh.faceIsCurved(faceChildren(:, 1)) = true;
            msh.elemIsCurved(elemChildren([1, 3])) = true;
        else
            msh.faceNormalx(faceChildren(:, 1)) = faceNormalxs(1);
            msh.faceNormaly(faceChildren(:, 1)) = faceNormalys(1);
            msh.faceSize(faceChildren(:, 1))    = faceSizes(1) / 2;
            msh.faceJac(faceChildren(:, 1))     = faceSizes(1) / 4;  
        end
        msh.faceType(faceChildren(:, 1))   = faceTypes(1);
        msh.faceLevel(faceChildren(:, 1))  = faceLevels(1) + 1;
        msh.faceParent(faceChildren(:, 1)) = elemFaces(1);      
        msh.nFaces = msh.nFaces + 2;
    else
        % It must be an internal or periodic face in this case
        nodeCoors(:, 1) = 0.5 * sum(elemNodeCoors(:, [4, 1]), 2);
        if (faceTypes(1) == 0)
            nodeIDs(1) = msh.faceNodes(2, faceChildren(1, 1));   
        elseif (faceTypes(1) == 1)
            if (faceNodes(1, 1) == elemNodes(4) && faceNodes(2, 1) == elemNodes(1))
                nodeIDs(1) = msh.faceNodes(2, faceChildren(1, 1));    
            elseif (faceNodes(1, 1) == elemNodes(1) && faceNodes(2, 1) == elemNodes(4))
                nodeIDs(1) = msh.faceNodes(2, faceChildren(1, 1));
            else            
                nodeIDs(1) = msh.nNodes + 1;
                msh.nodeCoor(:, nodeIDs(1)) = nodeCoors(:, 1);
                msh.nNodes = msh.nNodes + 1;
            end
        else
            error('Wrong face type')
        end   
        
        if (faceElems(1, 1) == ID)      
            subfIDs(:, 1) = faceChildren(:, 1);
            msh.faceElems(1, faceChildren(:, 1)) = elemChildren([3, 1]);
            msh.faceNums( 1, faceChildren(:, 1)) = 1; 
        else
            subfIDs(:, 1) = faceChildren([2, 1], 1);
            msh.faceElems(2, faceChildren(:, 1)) = elemChildren([1, 3]);
            msh.faceNums( 2, faceChildren(:, 1)) = 1; 
        end   
    end
    % second face
    if (faceChildren(1, 2) == 0)
        faceChildren(:, 2) = [msh.nFaces + 1; msh.nFaces + 2];
        msh.faceChildren(:, elemFaces(2)) = faceChildren(:, 2);
        
        nodeIDs(2) = msh.nNodes + 1;
        msh.nNodes = msh.nNodes + 1; 
        if (faceTypes(2) == 0)
            nodeCoors(:, 2) = 0.5 * sum(elemNodeCoors(:, 2 : 3), 2);
            nodeID = nodeIDs(2);
        elseif (faceTypes(2) == 1)
            nodeCoors(:, 2) = 0.5 * sum(elemNodeCoors(:, 2 : 3), 2);
            if (faceNodes(1, 2) == elemNodes(2) && faceNodes(2, 2) == elemNodes(3))
                nodeID = nodeIDs(2);
            elseif (faceNodes(1, 2) == elemNodes(3) && faceNodes(2, 2) == elemNodes(2))
                nodeID = nodeIDs(2);
            else
                nodeID = msh.nNodes + 1;
                msh.nodeCoor(:, nodeID) = 0.5 * sum(msh.nodeCoor(:, faceNodes(:, 2)), 2);
                msh.nNodes = msh.nNodes + 1; 
            end
        else
            nodeCoors(:, 2) = getBndNode(elemNodeCoors(:, 2), elemNodeCoors(:, 3), faceTypes(2));
            nodeID = nodeIDs(2);
        end
        msh.nodeCoor(:, nodeIDs(2)) = nodeCoors(:, 2);        

        msh.faceNodes(:, faceChildren(:, 2)) = [faceNodes(1, 2), nodeID; nodeID, faceNodes(2, 2)];        
        if (faceTypes(2) == 0 || faceTypes(2) == 1)                       
            if (faceElems(1, 2) == ID)               
                subfIDs(:, 2) = faceChildren(:, 2);
                msh.faceElems(1, faceChildren(:, 2)) = elemChildren([2, 4]);
                msh.faceElems(2, faceChildren(:, 2)) = faceElems(2, 2);
                msh.faceNums( 1, faceChildren(:, 2)) = 2;
                msh.faceNums( 2, faceChildren(:, 2)) = [2 * faceNums(2, 2) + 4, 2 * faceNums(2, 2) + 3];                                    
            else                
                subfIDs(:, 2) = faceChildren([2, 1], 2);
                msh.faceElems(1, faceChildren(:, 2)) = faceElems(1, 2);
                msh.faceElems(2, faceChildren(:, 2)) = elemChildren([4, 2]);
                msh.faceNums( 1, faceChildren(:, 2)) = [2 * faceNums(1, 2) + 3, 2 * faceNums(1, 2) + 4];
                msh.faceNums( 2, faceChildren(:, 2)) = 2;                              
            end               
        else
            subfIDs(:, 2) = faceChildren(:, 2);
            msh.faceElems(1, faceChildren(:, 2)) = elemChildren([2, 4]);
            msh.faceNums( 1, faceChildren(:, 2)) = 2;
        end   
        if faceIsCurved(2)
            tempNormals = [nodeCoors(2, 2) - elemNodeCoors(2, 2), elemNodeCoors(2, 3) - nodeCoors(2, 2); elemNodeCoors(1, 2) - nodeCoors(1, 2), nodeCoors(1, 2) - elemNodeCoors(1, 3)];
            tempSizes   = vecnorm(tempNormals);
            tempNormals = tempNormals ./ tempSizes;
            
            msh.faceNormalx(faceChildren(:, 2))  = tempNormals(1, :);
            msh.faceNormaly(faceChildren(:, 2))  = tempNormals(2, :);
            msh.faceSize(faceChildren(:, 2))     = tempSizes;
            msh.faceJac(faceChildren(:, 2))      = tempSizes / 2;  
            msh.faceIsCurved(faceChildren(:, 2)) = true;
            msh.elemIsCurved(elemChildren([2, 4])) = true;
        else
            msh.faceNormalx(faceChildren(:, 2)) = faceNormalxs(2);
            msh.faceNormaly(faceChildren(:, 2)) = faceNormalys(2);
            msh.faceSize(faceChildren(:, 2))    = faceSizes(2) / 2;
            msh.faceJac(faceChildren(:, 2))     = faceSizes(2) / 4;  
        end       
        msh.faceType(faceChildren(:, 2))    = faceTypes(2);
        msh.faceLevel(faceChildren(:, 2))   = faceLevels(2) + 1;
        msh.faceParent(faceChildren(:, 2))  = elemFaces(2);
        msh.nFaces = msh.nFaces + 2;        
    else
        % It must be an internal or periodic face in this case    
        nodeCoors(:, 2) = 0.5 * sum(elemNodeCoors(:, 2 : 3), 2);
        if (faceTypes(2) == 0)
            nodeIDs(2) = msh.faceNodes(2, faceChildren(1, 2));
        elseif (faceTypes(2) == 1)
            if (faceNodes(1, 2) == elemNodes(2) && faceNodes(2, 2) == elemNodes(3))
                nodeIDs(2) = msh.faceNodes(2, faceChildren(1, 2));
            elseif (faceNodes(1, 2) == elemNodes(3) && faceNodes(2, 2) == elemNodes(2))
                nodeIDs(2) = msh.faceNodes(2, faceChildren(1, 2));
            else            
                nodeIDs(2) = msh.nNodes + 1;
                msh.nodeCoor(:, nodeIDs(2)) = nodeCoors(:, 2);
                msh.nNodes = msh.nNodes + 1;  
            end
        else
            error('Wrong face type')
        end         
        if (faceElems(1, 2) == ID)
            subfIDs(:, 2) = faceChildren(:, 2);
            msh.faceElems(1, faceChildren(:, 2)) = elemChildren([2, 4]);
            msh.faceNums( 1, faceChildren(:, 2)) = 2;              
        else
            subfIDs(:, 2) = faceChildren([2, 1], 2);
            msh.faceElems(2, faceChildren(:, 2)) = elemChildren([4, 2]);
            msh.faceNums( 2, faceChildren(:, 2)) = 2;            
        end
    end
    % third face
    if (faceChildren(1, 3) == 0)
        faceChildren(:, 3) = [msh.nFaces + 1; msh.nFaces + 2];
        msh.faceChildren(:, elemFaces(3)) = faceChildren(:, 3);
        
        nodeIDs(3) = msh.nNodes + 1;
        msh.nNodes = msh.nNodes + 1;    
        if (faceTypes(3) == 0)
            nodeCoors(:, 3) = 0.5 * sum(elemNodeCoors(:, 1 : 2), 2); 
            nodeID = nodeIDs(3);
        elseif (faceTypes(3) == 1)
            nodeCoors(:, 3) = 0.5 * sum(elemNodeCoors(:, 1 : 2), 2); 
            if (faceNodes(1, 3) == elemNodes(1) && faceNodes(2, 3) == elemNodes(2))
                nodeID = nodeIDs(3);
            elseif (faceNodes(1, 3) == elemNodes(2) && faceNodes(2, 3) == elemNodes(1))
                nodeID = nodeIDs(3);
            else
                nodeID = msh.nNodes + 1;
                msh.nodeCoor(:, nodeID) = 0.5 * sum(msh.nodeCoor(:, faceNodes(:, 3)), 2);
                msh.nNodes = msh.nNodes + 1; 
            end
        else
            nodeCoors(:, 3) = getBndNode(elemNodeCoors(:, 1), elemNodeCoors(:, 2), faceTypes(3));
            nodeID = nodeIDs(3);
        end
        msh.nodeCoor(:, nodeIDs(3)) = nodeCoors(:, 3);         

        msh.faceNodes(:, faceChildren(:, 3)) = [faceNodes(1, 3), nodeID; nodeID, faceNodes(2, 3)];
        if (faceTypes(3) == 0 || faceTypes(3) == 1)                         
            if (faceElems(1, 3) == ID)                
                subfIDs(:, 3) = faceChildren(:, 3);
                msh.faceElems(1, faceChildren(:, 3)) = elemChildren([1, 2]);
                msh.faceElems(2, faceChildren(:, 3)) = faceElems(2, 3);
                msh.faceNums( 1, faceChildren(:, 3)) = 3;
                msh.faceNums( 2, faceChildren(:, 3)) = [2 * faceNums(2, 3) + 4, 2 * faceNums(2, 3) + 3];              
            else                
                subfIDs(:, 3) = faceChildren([2, 1], 3);
                msh.faceElems(1, faceChildren(:, 3)) = faceElems(1, 3);
                msh.faceElems(2, faceChildren(:, 3)) = elemChildren([2, 1]);
                msh.faceNums( 1, faceChildren(:, 3)) = [2 * faceNums(1, 3) + 3, 2 * faceNums(1, 3) + 4];
                msh.faceNums( 2, faceChildren(:, 3)) = 3;                   
            end               
        else
            subfIDs(:, 3) = faceChildren(:, 3);
            msh.faceElems(1, faceChildren(:, 3)) = elemChildren([1, 2]);
            msh.faceNums( 1, faceChildren(:, 3)) = 3;
        end  
        if faceIsCurved(3)
            tempNormals = [nodeCoors(2, 3) - elemNodeCoors(2, 1), elemNodeCoors(2, 2) - nodeCoors(2, 3); elemNodeCoors(1, 1) - nodeCoors(1, 3), nodeCoors(1, 3) - elemNodeCoors(1, 2)];
            tempSizes   = vecnorm(tempNormals);
            tempNormals = tempNormals ./ tempSizes;
            
            msh.faceNormalx(faceChildren(:, 3))  = tempNormals(1, :);
            msh.faceNormaly(faceChildren(:, 3))  = tempNormals(2, :);
            msh.faceSize(faceChildren(:, 3))     = tempSizes;
            msh.faceJac(faceChildren(:, 3))      = tempSizes / 2;      
            msh.faceIsCurved(faceChildren(:, 3)) = true;
            msh.elemIsCurved(elemChildren([1, 2])) = true;
        else
            msh.faceNormalx(faceChildren(:, 3)) = faceNormalxs(3);
            msh.faceNormaly(faceChildren(:, 3)) = faceNormalys(3);
            msh.faceSize(faceChildren(:, 3))    = faceSizes(3) / 2;
            msh.faceJac(faceChildren(:, 3))     = faceSizes(3) / 4;  
        end        
        msh.faceType(faceChildren(:, 3))    = faceTypes(3);
        msh.faceLevel(faceChildren(:, 3))   = faceLevels(3) + 1;
        msh.faceParent(faceChildren(:, 3))  = elemFaces(3);
        msh.nFaces = msh.nFaces + 2;        
    else
        % It must be an internal or periodic face in this case
        nodeCoors(:, 3) = 0.5 * sum(elemNodeCoors(:, 1 : 2), 2);
        if (faceTypes(3) == 0)
            nodeIDs(3) = msh.faceNodes(2, faceChildren(1, 3));
        elseif (faceTypes(3) == 1)
            if (faceNodes(1, 3) == elemNodes(1) && faceNodes(2, 3) == elemNodes(2))
                nodeIDs(3) = msh.faceNodes(2, faceChildren(1, 3));
            elseif (faceNodes(1, 3) == elemNodes(2) && faceNodes(2, 3) == elemNodes(1))
                nodeIDs(3) = msh.faceNodes(2, faceChildren(1, 3));
            else            
                nodeIDs(3) = msh.nNodes + 1;
                msh.nodeCoor(:, nodeIDs(3)) = nodeCoors(:, 3);
                msh.nNodes = msh.nNodes + 1; 
            end
        else
            error('Wrong face type')
        end         
        if (faceElems(1, 3) == ID)           
            subfIDs(:, 3) = faceChildren(:, 3);
            msh.faceElems(1, faceChildren(:, 3)) = elemChildren([1, 2]);
            msh.faceNums( 1, faceChildren(:, 3)) = 3;               
        else
            subfIDs(:, 3) = faceChildren([2, 1], 3);
            msh.faceElems(2, faceChildren(:, 3)) = elemChildren([2, 1]);
            msh.faceNums( 2, faceChildren(:, 3)) = 3;             
        end         
    end   
    % fourth face
    if (faceChildren(1, 4) == 0)
        faceChildren(:, 4) = [msh.nFaces + 1; msh.nFaces + 2];
        msh.faceChildren(:, elemFaces(4)) = faceChildren(:, 4);
        
        nodeIDs(4) = msh.nNodes + 1;
        msh.nNodes = msh.nNodes + 1; 
        if (faceTypes(4) == 0)
            nodeCoors(:, 4) = 0.5 * sum(elemNodeCoors(:, 3 : 4), 2); 
            nodeID = nodeIDs(4);
        elseif (faceTypes(4) == 1)
            nodeCoors(:, 4) = 0.5 * sum(elemNodeCoors(:, 3 : 4), 2); 
            if (faceNodes(1, 4) == elemNodes(3) && faceNodes(2, 4) == elemNodes(4))
                nodeID = nodeIDs(4);
            elseif (faceNodes(1, 4) == elemNodes(4) && faceNodes(2, 4) == elemNodes(3))
                nodeID = nodeIDs(4);
            else
                nodeID = msh.nNodes + 1;
                msh.nodeCoor(:, nodeID) = 0.5 * sum(msh.nodeCoor(:, faceNodes(:, 4)), 2);
                msh.nNodes = msh.nNodes + 1; 
            end
        else
            nodeCoors(:, 4) = getBndNode(elemNodeCoors(:, 3), elemNodeCoors(:, 4), faceTypes(4));
            nodeID = nodeIDs(4);
        end
        msh.nodeCoor(:, nodeIDs(4)) = nodeCoors(:, 4);         

        msh.faceNodes(:, faceChildren(:, 4)) = [faceNodes(1, 4), nodeID; nodeID, faceNodes(2, 4)];
        if (faceTypes(4) == 0 || faceTypes(4) == 1)                    
            if (faceElems(1, 4) == ID)             
                subfIDs(:, 4) = faceChildren(:, 4);
                msh.faceElems(1, faceChildren(:, 4)) = elemChildren([4, 3]);
                msh.faceElems(2, faceChildren(:, 4)) = faceElems(2, 4);
                msh.faceNums( 1, faceChildren(:, 4)) = 4;
                msh.faceNums( 2, faceChildren(:, 4)) = [2 * faceNums(2, 4) + 4, 2 * faceNums(2, 4) + 3];                
            else   
                subfIDs(:, 4) = faceChildren([2, 1], 4);
                msh.faceElems(1, faceChildren(:, 4)) = faceElems(1, 4);
                msh.faceElems(2, faceChildren(:, 4)) = elemChildren([3, 4]);
                msh.faceNums( 1, faceChildren(:, 4)) = [2 * faceNums(1, 4) + 3, 2 * faceNums(1, 4) + 4];
                msh.faceNums( 2, faceChildren(:, 4)) = 4;                 
            end               
        else
            subfIDs(:, 4) = faceChildren(:, 4);
            msh.faceElems(1, faceChildren(:, 4)) = elemChildren([4, 3]);
            msh.faceNums( 1, faceChildren(:, 4)) = 4;            
        end
        if faceIsCurved(4)
            tempNormals = [nodeCoors(2, 4) - elemNodeCoors(2, 3), elemNodeCoors(2, 4) - nodeCoors(2, 4); elemNodeCoors(1, 3) - nodeCoors(1, 4), nodeCoors(1, 4) - elemNodeCoors(1, 4)];
            tempSizes   = vecnorm(tempNormals);
            tempNormals = tempNormals ./ tempSizes;
            
            msh.faceNormalx(faceChildren(:, 4))  = tempNormals(1, :);
            msh.faceNormaly(faceChildren(:, 4))  = tempNormals(2, :);
            msh.faceSize(faceChildren(:, 4))     = tempSizes;
            msh.faceJac(faceChildren(:, 4))      = tempSizes / 2;
            msh.faceIsCurved(faceChildren(:, 4)) = true;
            msh.elemIsCurved(elemChildren([3, 4])) = true;
        else
            msh.faceNormalx(faceChildren(:, 4)) = faceNormalxs(4);
            msh.faceNormaly(faceChildren(:, 4)) = faceNormalys(4);
            msh.faceSize(faceChildren(:, 4))    = faceSizes(4) / 2;
            msh.faceJac(faceChildren(:, 4))     = faceSizes(4) / 4;  
        end                
        msh.faceType(faceChildren(:, 4))    = faceTypes(4);
        msh.faceLevel(faceChildren(:, 4))   = faceLevels(4) + 1;
        msh.faceParent(faceChildren(:, 4))  = elemFaces(4);
        msh.nFaces = msh.nFaces + 2;        
    else
        % It must be an internal or periodic face in this case
        nodeCoors(:, 4) = 0.5 * sum(elemNodeCoors(:, 3 : 4), 2);
        if (faceTypes(4) == 0)
            nodeIDs(4) = msh.faceNodes(2, faceChildren(1, 4));
        elseif (faceTypes(4) == 1)
            if (faceNodes(1, 4) == elemNodes(3) && faceNodes(2, 4) == elemNodes(4))
                nodeIDs(4) = msh.faceNodes(2, faceChildren(1, 4));
            elseif (faceNodes(1, 4) == elemNodes(4) && faceNodes(2, 4) == elemNodes(3))
                nodeIDs(4) = msh.faceNodes(2, faceChildren(1, 4));
            else            
                nodeIDs(4) = msh.nNodes + 1;
                msh.nodeCoor(:, nodeIDs(4)) = nodeCoors(:, 4);
                msh.nNodes = msh.nNodes + 1;   
            end
        else
            error('Wrong face type')
        end         
        if (faceElems(1, 4) == ID)
            subfIDs(:, 4) = faceChildren(:, 4);
            msh.faceElems(1, faceChildren(:, 4)) = elemChildren([4, 3]);
            msh.faceNums( 1, faceChildren(:, 4)) = 4;             
        else
            subfIDs(:, 4) = faceChildren([2, 1], 4);
            msh.faceElems(2, faceChildren(:, 4)) = elemChildren([3, 4]);
            msh.faceNums( 2, faceChildren(:, 4)) = 4;              
        end         
    end   
    % interior nodes and four interior faces
    nodeCoors(:, 5) = 0.25 * sum(nodeCoors(:, 1 : 4), 2);
    nodeIDs(5) = msh.nNodes + 1;
    msh.nodeCoor(:, nodeIDs(5)) = nodeCoors(:, 5);
    msh.nNodes = msh.nNodes + 1;
    
    subfIDs(:, 5) = [msh.nFaces + 1; msh.nFaces + 2]; 
    msh.faceNodes(:, subfIDs(:, 5)) = [nodeIDs(3), nodeIDs(5); nodeIDs(5), nodeIDs(4)];      
    msh.faceElems(1, subfIDs(:, 5)) = elemChildren([1, 3]);
    msh.faceElems(2, subfIDs(:, 5)) = elemChildren([2, 4]);
    msh.faceNums( 1, subfIDs(:, 5)) = 2;
    msh.faceNums( 2, subfIDs(:, 5)) = 1;  
    tempNormals = [nodeCoors(2, 5) - nodeCoors(2, 3), nodeCoors(2, 4) - nodeCoors(2, 5); nodeCoors(1, 3) - nodeCoors(1, 5), nodeCoors(1, 5) - nodeCoors(1, 4)];
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    msh.faceNormalx(subfIDs(:, 5)) = tempNormals(1, :);
    msh.faceNormaly(subfIDs(:, 5)) = tempNormals(2, :);
    msh.faceSize(subfIDs(:, 5))    = tempSizes;
    msh.faceJac(subfIDs(:, 5))     = tempSizes / 2;  
    subfIDs(:, 6) = [msh.nFaces + 3; msh.nFaces + 4]; 
    msh.faceNodes(:, subfIDs(:, 6)) = [nodeIDs(2), nodeIDs(5); nodeIDs(5), nodeIDs(1)];
    msh.faceElems(1, subfIDs(:, 6)) = elemChildren([2, 1]);
    msh.faceElems(2, subfIDs(:, 6)) = elemChildren([4, 3]);
    msh.faceNums( 1, subfIDs(:, 6)) = 4;
    msh.faceNums( 2, subfIDs(:, 6)) = 3;  
    tempNormals = [nodeCoors(2, 5) - nodeCoors(2, 2), nodeCoors(2, 1) - nodeCoors(2, 5); nodeCoors(1, 2) - nodeCoors(1, 5), nodeCoors(1, 5) - nodeCoors(1, 1)];
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    msh.faceNormalx(subfIDs(:, 6)) = tempNormals(1, :);
    msh.faceNormaly(subfIDs(:, 6)) = tempNormals(2, :);
    msh.faceSize(subfIDs(:, 6))    = tempSizes;
    msh.faceJac(subfIDs(:, 6))     = tempSizes / 2;     
    msh.nFaces = msh.nFaces + 4;
    
    % element information     
    subElemNodeCoorx = [elemNodeCoors(1, 1), nodeCoors(1, 3),     nodeCoors(1, 1),     nodeCoors(1, 5);...
                        nodeCoors(1, 3),     elemNodeCoors(1, 2), nodeCoors(1, 5),     nodeCoors(1, 2);...
                        nodeCoors(1, 5),     nodeCoors(1, 2),     nodeCoors(1, 4),     elemNodeCoors(1, 3);...
                        nodeCoors(1, 1),     nodeCoors(1, 5),     elemNodeCoors(1, 4), nodeCoors(1, 4)]; 
                    
    subElemNodeCoory = [elemNodeCoors(2, 1), nodeCoors(2, 3),     nodeCoors(2, 1),     nodeCoors(2, 5);...
                        nodeCoors(2, 3),     elemNodeCoors(2, 2), nodeCoors(2, 5),     nodeCoors(2, 2);...
                        nodeCoors(2, 5),     nodeCoors(2, 2),     nodeCoors(2, 4),     elemNodeCoors(2, 3);...
                        nodeCoors(2, 1),     nodeCoors(2, 5),     elemNodeCoors(2, 4), nodeCoors(2, 4)];   
                    
    x    = quad.N     * subElemNodeCoorx;
    xxi  = quad.N_xi  * subElemNodeCoorx;
    xeta = quad.N_eta * subElemNodeCoorx;
    y    = quad.N     * subElemNodeCoory;
    yxi  = quad.N_xi  * subElemNodeCoory;
    yeta = quad.N_eta * subElemNodeCoory;     
    J    = xxi .* yeta - xeta .* yxi;
                    
    si = 0.5 * ((subElemNodeCoorx(1, :) - subElemNodeCoorx(3, :)) .* (subElemNodeCoory(2, :) - subElemNodeCoory(4, :)) - (subElemNodeCoorx(2, :) - subElemNodeCoorx(4, :)) .* (subElemNodeCoory(1, :) - subElemNodeCoory(3, :)));
    msh.elemSize(elemChildren)     = si;
    msh.elemDiam(elemChildren)     = sqrt(si);
    msh.elemNodes(:, elemChildren) = [elemNodes(1), nodeIDs(3),   nodeIDs(1),   nodeIDs(5);...
                                      nodeIDs(3),   elemNodes(2), nodeIDs(5),   nodeIDs(2);...
                                      nodeIDs(5),   nodeIDs(2),   nodeIDs(4),   elemNodes(3);...
                                      nodeIDs(1),   nodeIDs(5),   elemNodes(4), nodeIDs(4)];
    msh.elemFaces(:, elemChildren) = [subfIDs(2, 1), subfIDs(1, 5), subfIDs(1, 1), subfIDs(2, 5);...
                                      subfIDs(1, 5), subfIDs(1, 2), subfIDs(2, 5), subfIDs(2, 2);...
                                      subfIDs(1, 3), subfIDs(2, 3), subfIDs(2, 6), subfIDs(1, 6);...
                                      subfIDs(2, 6), subfIDs(1, 6), subfIDs(2, 4), subfIDs(1, 4)]; 
    msh.elemGPx(:, elemChildren)   = x;
    msh.elemGPy(:, elemChildren)   = y;
    msh.elemJac(:, elemChildren)   = J;
    msh.elemJxix(:, elemChildren)  = yeta; 
    msh.elemJxiy(:, elemChildren)  = -xeta; 
    msh.elemJetax(:, elemChildren) = -yxi; 
    msh.elemJetay(:, elemChildren) = xxi;         
    msh.nElems = msh.nElems + 4;

    % inverse of mass matrix
    mi1 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 1)));
    mi2 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 2)));
    mi3 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 3)));
    mi4 = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 4)));
    msh.massMatInv{elemChildren(1)} = mi1;
    msh.massMatInv{elemChildren(2)} = mi2;
    msh.massMatInv{elemChildren(3)} = mi3;
    msh.massMatInv{elemChildren(4)} = mi4;
    
    % Project the solution in the coarse element to the four refined
    % elements
    UE = UT(:, nv * (ID - 1) + 1 : nv * ID);
    if msh.elemIsCurved(ID) 
        [xi, eta] = phys2Ref_quad(x, y, repmat(elemNodeCoors(1, :)', [1, 4]), repmat(elemNodeCoors(2, :)', [1, 4]));
        UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) = mi1 * (bs.phitw{1} * (basisFunctionSet_square(xi(:, 1), eta(:, 1), bs.type, [0, 0]) * UE .* J(:, 1)));   
        UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) = mi2 * (bs.phitw{1} * (basisFunctionSet_square(xi(:, 2), eta(:, 2), bs.type, [0, 0]) * UE .* J(:, 2)));
        UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) = mi3 * (bs.phitw{1} * (basisFunctionSet_square(xi(:, 3), eta(:, 3), bs.type, [0, 0]) * UE .* J(:, 3)));
        UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4)) = mi4 * (bs.phitw{1} * (basisFunctionSet_square(xi(:, 4), eta(:, 4), bs.type, [0, 0]) * UE .* J(:, 4)));        
    else
        UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(1)) = mi1 * (bs.phitw{1} * (bs.phi_sube{1} * UE .* J(:, 1)));   
        UT(:, nv * (elemChildren(2) - 1) + 1 : nv * elemChildren(2)) = mi2 * (bs.phitw{1} * (bs.phi_sube{2} * UE .* J(:, 2)));
        UT(:, nv * (elemChildren(3) - 1) + 1 : nv * elemChildren(3)) = mi3 * (bs.phitw{1} * (bs.phi_sube{3} * UE .* J(:, 3)));
        UT(:, nv * (elemChildren(4) - 1) + 1 : nv * elemChildren(4)) = mi4 * (bs.phitw{1} * (bs.phi_sube{4} * UE .* J(:, 4)));
    end
end    
msh.faceNormalx(:, msh.nFaces + 1 : end)  = []; 
msh.faceNormaly(:, msh.nFaces + 1 : end)  = []; 
msh.faceSize(msh.nFaces + 1 : end)        = []; 
msh.faceType(msh.nFaces + 1 : end)        = [];
msh.faceNodes(:, msh.nFaces + 1 : end)    = [];
msh.faceElems(:, msh.nFaces + 1 : end)    = []; 
msh.faceNums(:, msh.nFaces + 1 : end)     = []; 
msh.faceLevel(msh.nFaces + 1 : end)       = []; 
msh.faceParent(msh.nFaces + 1 : end)      = []; 
msh.faceChildren(:, msh.nFaces + 1 : end) = []; 
msh.faceJac(:, msh.nFaces + 1 : end)      = []; 
msh.faceIsCurved(msh.nFaces + 1 : end)    = []; 
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

% Make sure all the four siblings need to be coarsened and store the IDs of
% their parents
npe = 0;
peIDs = zeros(1, floor(nce / 4));
i = 1;
while i <= nce
    ID = ceIDs(i);
    parentID = msh.elemParent(ID);
    % Make sure all the four siblings need to be coarsened
    siblingIDs = msh.elemChildren(:, parentID);  
    if any(~shouldCoarsen(siblingIDs))
        shouldCoarsen(ID) = false;
        i = i + 1;
        continue;
    end
    % Store the parent ID
    npe = npe + 1;
    peIDs(npe) = parentID;
    
    i = i + 4;
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

% Update face-element connections, or delete elements, faces, nodes and 
% store their IDs 
deIDs = zeros(1, 4 * npe);
dfIDs = zeros(1, 12 * npe);
dnIDs = zeros(1, 5 * npe);
nde = 0;
ndf = 0;
ndn = 0;
for i = 1 : npe
    ID = peIDs(i);
    elemNodes    = msh.elemNodes(:, ID);
    elemFaces    = msh.elemFaces(:, ID);
    elemChildren = msh.elemChildren(:, ID);  
    faceTypes    = msh.faceType(elemFaces);
    faceNodes    = msh.faceNodes(:, elemFaces);
    faceElems    = msh.faceElems(:, elemFaces);
    faceNums     = msh.faceNums(:, elemFaces);
    faceChildren = msh.faceChildren(:, elemFaces);
    
    % Make sure the level difference of adjacent elements is less than 1
    if any(msh.faceChildren(1, faceChildren) ~= 0)
        shouldCoarsen(elemChildren) = false;
        continue;  
    end
    
    % Project the solution in the four refined elements into the coarse
    % element
    J = msh.elemJac(:, elemChildren);
    u4 = bs.phi{1} * UT(:, nv * (elemChildren(1) - 1) + 1 : nv * elemChildren(4));
    if msh.elemIsCurved(ID)    
        elemNodeCoors = msh.nodeCoor(:, elemNodes);
        x = msh.elemGPx(:, elemChildren);
        y = msh.elemGPy(:, elemChildren);
        [xi, eta] = phys2Ref_quad(x, y, repmat(elemNodeCoors(1, :)', [1, 4]), repmat(elemNodeCoors(2, :)', [1, 4]));
        
        UT(:, nv * (ID - 1) + 1 : nv * ID) = msh.massMatInv{ID} * ((quad.weights .* basisFunctionSet_square(xi(:, 1), eta(:, 1), bs.type, [0, 0]))' * (u4(:, 1 : nv) .* J(:, 1)) ...
                                                                 + (quad.weights .* basisFunctionSet_square(xi(:, 2), eta(:, 2), bs.type, [0, 0]))' * (u4(:, nv + 1 : 2 * nv) .* J(:, 2)) ...
                                                                 + (quad.weights .* basisFunctionSet_square(xi(:, 3), eta(:, 3), bs.type, [0, 0]))' * (u4(:, 2 * nv + 1 : 3 * nv) .* J(:, 3)) ...
                                                                 + (quad.weights .* basisFunctionSet_square(xi(:, 4), eta(:, 4), bs.type, [0, 0]))' * (u4(:, 3 * nv + 1 : end) .* J(:, 4)));               
    else
       UT(:, nv * (ID - 1) + 1 : nv * ID) = msh.massMatInv{ID} * (bs.phitw_sube{1} * (u4(:, 1 : nv) .* J(:, 1)) ...
                                                                + bs.phitw_sube{2} * (u4(:, nv + 1 : 2 * nv) .* J(:, 2)) ...
                                                                + bs.phitw_sube{3} * (u4(:, 2 * nv + 1 : 3 * nv) .* J(:, 3)) ...
                                                                + bs.phitw_sube{4} * (u4(:, 3 * nv + 1 : end) .* J(:, 4)));   
    end
                                           
    % Update face-element connections or delete subfaces and nodes
    % first face
    if (faceTypes(1) == 0 || faceTypes(1) == 1)
        isLeft = faceElems(1, 1) == ID;
        ID1 = faceElems(1, 1) * ~isLeft + faceElems(2, 1) * isLeft;
        if (msh.elemChildren(1, ID1) ~= 0)
            if isLeft
                msh.faceElems(1, faceChildren(:, 1)) = ID;
                msh.faceNums( 1, faceChildren(:, 1)) = [2 * faceNums(1, 1) + 3, 2 * faceNums(1, 1) + 4];
                if ~(faceNodes(1, 1) == elemNodes(4) && faceNodes(2, 1) == elemNodes(1))
                    dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(1));
                    ndn = ndn + 1;
                end
            else
                msh.faceElems(2, faceChildren(:, 1)) = ID;
                msh.faceNums( 2, faceChildren(:, 1)) = [2 * faceNums(2, 1) + 4, 2 * faceNums(2, 1) + 3];
                if ~(faceNodes(1, 1) == elemNodes(1) && faceNodes(2, 1) == elemNodes(4))
                    dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(1));
                    ndn = ndn + 1;
                end                
            end
        else
            msh.faceChildren(:, elemFaces(1)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 1)';
            ndf = ndf + 2;

            dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(1));
            ndn = ndn + 1;
            if isLeft
                if ~(faceNodes(1, 1) == elemNodes(4) && faceNodes(2, 1) == elemNodes(1))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 1));
                    ndn = ndn + 1;
                end
            else
                if ~(faceNodes(1, 1) == elemNodes(1) && faceNodes(2, 1) == elemNodes(4))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 1));
                    ndn = ndn + 1;
                end                
            end
        end
    else        
        msh.faceChildren(:, elemFaces(1)) = 0;
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 1)';
        ndf = ndf + 2;

        dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(1));
        ndn = ndn + 1;
    end
    % second face
    if (faceTypes(2) == 0 || faceTypes(2) == 1)
        isLeft = faceElems(1, 2) == ID;
        ID2 = faceElems(1, 2) * ~isLeft + faceElems(2, 2) * isLeft;
        if (msh.elemChildren(1, ID2) ~= 0)
            if isLeft
                msh.faceElems(1, faceChildren(:, 2)) = ID;
                msh.faceNums( 1, faceChildren(:, 2)) = [2 * faceNums(1, 2) + 3, 2 * faceNums(1, 2) + 4];
                if ~(faceNodes(1, 2) == elemNodes(2) && faceNodes(2, 2) == elemNodes(3))
                    dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(4));
                    ndn = ndn + 1;
                end                
            else
                msh.faceElems(2, faceChildren(:, 2)) = ID;
                msh.faceNums( 2, faceChildren(:, 2)) = [2 * faceNums(2, 2) + 4, 2 * faceNums(2, 2) + 3]; 
                if ~(faceNodes(1, 2) == elemNodes(3) && faceNodes(2, 2) == elemNodes(2))
                    dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(4));
                    ndn = ndn + 1;
                end                 
            end
        else
            msh.faceChildren(:, elemFaces(2)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 2)';
            ndf = ndf + 2;

            dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(4));
            ndn = ndn + 1; 
            if isLeft
                if ~(faceNodes(1, 2) == elemNodes(2) && faceNodes(2, 2) == elemNodes(3))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 2));
                    ndn = ndn + 1;
                end
            else
                if ~(faceNodes(1, 2) == elemNodes(3) && faceNodes(2, 2) == elemNodes(2))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 2));
                    ndn = ndn + 1;                    
                end
            end
        end
    else
        msh.faceChildren(:, elemFaces(2)) = 0;     
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 2)';
        ndf = ndf + 2;

        dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(4));
        ndn = ndn + 1;     
    end   
    % third face
    if (faceTypes(3) == 0 || faceTypes(3) == 1)
        isLeft = faceElems(1, 3) == ID;
        ID3 = faceElems(1, 3) * ~isLeft + faceElems(2, 3) * isLeft;
        if (msh.elemChildren(1, ID3) ~= 0)
            if isLeft
                msh.faceElems(1, faceChildren(:, 3)) = ID;
                msh.faceNums( 1, faceChildren(:, 3)) = [2 * faceNums(1, 3) + 3, 2 * faceNums(1, 3) + 4];
                if ~(faceNodes(1, 3) == elemNodes(1) && faceNodes(2, 3) == elemNodes(2))
                    dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(1));
                    ndn = ndn + 1; 
                end
            else
                msh.faceElems(2, faceChildren(:, 3)) = ID;
                msh.faceNums( 2, faceChildren(:, 3)) = [2 * faceNums(2, 3) + 4, 2 * faceNums(2, 3) + 3];
                if ~(faceNodes(1, 3) == elemNodes(2) && faceNodes(2, 3) == elemNodes(1))
                    dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(1));
                    ndn = ndn + 1; 
                end                
            end
        else
            msh.faceChildren(:, elemFaces(3)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 3)';
            ndf = ndf + 2;

            dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(1));
            ndn = ndn + 1;   
            if isLeft
                if ~(faceNodes(1, 3) == elemNodes(1) && faceNodes(2, 3) == elemNodes(2))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 3));
                    ndn = ndn + 1;
                end
            else
                if ~(faceNodes(1, 3) == elemNodes(2) && faceNodes(2, 3) == elemNodes(1))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 3));
                    ndn = ndn + 1;                    
                end
            end
        end
    else
        msh.faceChildren(:, elemFaces(3)) = 0;    
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 3)';
        ndf = ndf + 2;

        dnIDs(ndn + 1) = msh.elemNodes(2, elemChildren(1));
        ndn = ndn + 1;
    end    
    % fourth face
    if (faceTypes(4) == 0 || faceTypes(4) == 1)
        isLeft = faceElems(1, 4) == ID;
        ID4 = faceElems(1, 4) * ~isLeft + faceElems(2, 4) * isLeft;
        if (msh.elemChildren(1, ID4) ~= 0)
            if isLeft
                msh.faceElems(1, faceChildren(:, 4)) = ID;
                msh.faceNums( 1, faceChildren(:, 4)) = [2 * faceNums(1, 4) + 3, 2 * faceNums(1, 4) + 4];
                if ~(faceNodes(1, 4) == elemNodes(3) && faceNodes(2, 4) == elemNodes(4))
                    dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(4));
                    ndn = ndn + 1;  
                end
            else
                msh.faceElems(2, faceChildren(:, 4)) = ID;
                msh.faceNums( 2, faceChildren(:, 4)) = [2 * faceNums(2, 4) + 4, 2 * faceNums(2, 4) + 3];
                if ~(faceNodes(1, 4) == elemNodes(4) && faceNodes(2, 4) == elemNodes(3))
                    dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(4));
                    ndn = ndn + 1;  
                end                
            end
        else          
            msh.faceChildren(:, elemFaces(4)) = 0;
            dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 4)';
            ndf = ndf + 2;

            dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(4));
            ndn = ndn + 1;     
            if isLeft
                if ~(faceNodes(1, 4) == elemNodes(3) && faceNodes(2, 4) == elemNodes(4))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 4));
                    ndn = ndn + 1;                    
                end
            else
                if ~(faceNodes(1, 4) == elemNodes(4) && faceNodes(2, 4) == elemNodes(3))
                    dnIDs(ndn + 1) = msh.faceNodes(2, faceChildren(1, 4));
                    ndn = ndn + 1;                      
                end
            end
        end
    else
        msh.faceChildren(:, elemFaces(4)) = 0;     
        dfIDs(ndf + 1 : ndf + 2) = faceChildren(:, 4)';
        ndf = ndf + 2;

        dnIDs(ndn + 1) = msh.elemNodes(4, elemChildren(4));
        ndn = ndn + 1;   
    end  
    % interior node and four interior faces
    childFaceID = msh.elemFaces(2, elemChildren(1));
    dfIDs(ndf + 1 : ndf + 4) = childFaceID : childFaceID + 3;
    ndf = ndf + 4;
    
    dnIDs(ndn + 1) = msh.elemNodes(3, elemChildren(1));
    ndn = ndn + 1;
    
    % Delete elements
    msh.elemChildren(:, ID) = 0;
    deIDs(nde + 1 : nde + 4) = elemChildren;
    nde = nde + 4;
end
deIDs(nde + 1 : end) = [];
dfIDs(ndf + 1 : end) = [];
dnIDs(ndn + 1 : end) = [];

% Re-index the remaining elements, faces and nodes
remEIDs = 1 : msh.nElems;
remFIDs = 1 : msh.nFaces;
remNIDs = 1 : msh.nNodes;
remEIDs(deIDs) = []; % remaining element old IDs
remFIDs(dfIDs) = []; % remaining face old IDs
remNIDs(dnIDs) = []; % remaining node old IDs

neIDs = zeros(1, msh.nElems); 
nfIDs = zeros(1, msh.nFaces); 
nnIDs = zeros(1, msh.nNodes);
neIDs(remEIDs) = 1 : length(remEIDs); % new element IDs
nfIDs(remFIDs) = 1 : length(remFIDs); % new face IDs
nnIDs(remNIDs) = 1 : length(remNIDs); % new node IDs

% Refresh element-node, element-face connections, parent and child IDs of elements
msh.elemNodes(:, remEIDs) = nnIDs(msh.elemNodes(:, remEIDs));
msh.elemFaces(:, remEIDs) = nfIDs(msh.elemFaces(:, remEIDs));

remEIDs_hasParent = remEIDs(msh.elemParent(remEIDs) ~= 0);
msh.elemParent(remEIDs_hasParent) = neIDs(msh.elemParent(remEIDs_hasParent));

remEIDs_hasChild = remEIDs(msh.elemChildren(1, remEIDs) ~= 0);
msh.elemChildren(:, remEIDs_hasChild) = neIDs(msh.elemChildren(:, remEIDs_hasChild));

% Refresh face-element, face-node connections, parent and child IDs of faces
msh.faceNodes(:, remFIDs) = nnIDs(msh.faceNodes(:, remFIDs));

msh.faceElems(1, remFIDs) = neIDs(msh.faceElems(1, remFIDs));
remFIDs_int = remFIDs(msh.faceType(remFIDs) == 0 | msh.faceType(remFIDs) == 1);
msh.faceElems(2, remFIDs_int) = neIDs(msh.faceElems(2, remFIDs_int));

remFIDs_hasParent = remFIDs(msh.faceParent(remFIDs) ~= 0);
msh.faceParent(remFIDs_hasParent) = nfIDs(msh.faceParent(remFIDs_hasParent));

remFIDs_hasChild = remFIDs(msh.faceChildren(1, remFIDs) ~= 0);
msh.faceChildren(:, remFIDs_hasChild) = nfIDs(msh.faceChildren(:, remFIDs_hasChild));

% Finally remove the deleted elements and faces from the mesh
msh.elemSize(deIDs)        = [];
msh.elemDiam(deIDs)        = [];
msh.elemNodes(:, deIDs)    = [];
msh.elemFaces(:, deIDs)    = [];
msh.elemLevel(deIDs)       = [];
msh.elemParent(deIDs)      = [];
msh.elemChildren(:, deIDs) = [];
msh.elemGPx(:, deIDs)      = [];
msh.elemGPy(:, deIDs)      = [];
msh.elemJac(:, deIDs)      = [];
msh.elemJxix(:, deIDs)     = [];
msh.elemJxiy(:, deIDs)     = [];
msh.elemJetax(:, deIDs)    = [];
msh.elemJetay(:, deIDs)    = [];
msh.elemIsCurved(deIDs)    = [];
msh.faceNormalx(:, dfIDs)  = [];
msh.faceNormaly(:, dfIDs)  = [];
msh.faceSize(dfIDs)        = [];
msh.faceType(dfIDs)        = [];
msh.faceNodes(:, dfIDs)    = [];
msh.faceElems(:, dfIDs)    = [];
msh.faceNums(:, dfIDs)     = [];
msh.faceLevel(dfIDs)       = []; 
msh.faceParent(dfIDs)      = []; 
msh.faceChildren(:, dfIDs) = []; 
msh.faceJac(:, dfIDs)      = [];
msh.faceIsCurved(dfIDs)    = [];
msh.nodeCoor(:, dnIDs)     = [];
msh.massMatInv(deIDs)      = [];
msh.nElems = length(remEIDs);
msh.nFaces = length(remFIDs);
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

end




