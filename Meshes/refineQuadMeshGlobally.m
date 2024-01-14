% Copyright (C) Xiangyi Meng.

function msh = refineQuadMeshGlobally(msh0, getBndNode, quad, bs)

if (msh0.type ~= 211) && (msh0.type ~= 212)
    error('Wrong mesh type')
end
if (msh0.nElems ~= msh0.nLElems)
    error('All the elements of the background mesh should be leaf elements')
end

if (nargin < 2) || isempty(getBndNode)
    getBndNode = @(p1, p2, faceType)0.5 * (p1 + p2);
end

if (nargin < 3) || isempty(quad)
    quad = GaussQuadratureRule_quad([2, 2], 201);
end

if (nargin < 4) || isempty(bs)
    quad1 = GaussQuadratureRule_line(2, 101);
    bs = setBasisFunctionSet_square(quad1, quad, 201, 2);
end
%**************************************************************************
msh          = getEmptyMesh;
msh.dm       = msh0.dm;
msh.N        = 2 * msh0.N;
msh.type     = msh0.type;
msh.maxLevel = msh0.maxLevel;
msh.bndTypes = msh0.bndTypes;
msh.nElems   = 0;
msh.nFaces   = 0;
msh.nNodes   = msh0.nNodes;
msh.nLElems  = 4 * msh0.nElems; 
msh.nLFaces  = 4 * msh0.nElems + 2 * msh0.nFaces; 
msh.LElems   = 1 : msh.nLElems;
msh.LFaces   = 1 : msh.nLFaces;
msh.nGPs     = msh0.nGPs;
%**************************************************************************
si = 4 * msh0.nElems;
msh.elemSize     = zeros(1, si);
msh.elemDiam     = zeros(1, si);
msh.elemNodes    = zeros(4, si);
msh.elemFaces    = zeros(4, si);
msh.elemGPx      = zeros(quad.np, si);
msh.elemGPy      = zeros(quad.np, si);
msh.elemJac      = zeros(quad.np, si);
msh.elemJxix     = zeros(quad.np, si); 
msh.elemJxiy     = zeros(quad.np, si); 
msh.elemJetax    = zeros(quad.np, si); 
msh.elemJetay    = zeros(quad.np, si); 
msh.elemIsCurved = false(1, si);
si = 4 * msh0.nElems + 2 * msh0.nFaces;
msh.faceNormalx  = zeros(1, si); 
msh.faceNormaly  = zeros(1, si);
msh.faceSize     = zeros(1, si);
msh.faceType     = zeros(1, si);
msh.faceNodes    = zeros(2, si);
msh.faceElems    = zeros(2, si); 
msh.faceNums     = zeros(2, si); 
msh.faceJac      = zeros(1, si); 
msh.faceIsCurved = false(1, si); 
si = msh0.nNodes + msh0.nElems + msh0.nFaces + sum(msh0.N);
msh.nodeCoor = zeros(2, si);
msh.nodeCoor(:, 1 : msh0.nNodes) = msh0.nodeCoor;
%**************************************************************************
hasVisited   = false(1, msh0.nFaces);
totalFaceIDs = zeros(2, msh0.nFaces);
newFaceIDs   = zeros(2, 6);
newNodeIDs   = zeros(1, 5);
newNodeCoors = zeros(2, 5);
for ID = 1 : msh0.nElems
    % current loop element information
    elemNodes     = msh0.elemNodes(:, ID);
    elemFaces     = msh0.elemFaces(:, ID);
    elemNodeCoors = msh0.nodeCoor(:, elemNodes);
    faceNormalxs  = msh0.faceNormalx(elemFaces);
    faceNormalys  = msh0.faceNormaly(elemFaces);
    faceSizes     = msh0.faceSize(elemFaces);
    faceTypes     = msh0.faceType(elemFaces);
    faceNodes     = msh0.faceNodes(:, elemFaces);
    faceElems     = msh0.faceElems(:, elemFaces);
    faceNums      = msh0.faceNums(:, elemFaces);
    faceIsCurved  = msh0.faceIsCurved(elemFaces);
    
    % creation of nodes, faces and elements of the new mesh
    for i = 1 : 4
        localFaceNodeIDs = getLocalFaceNodeIDs(i); 
        faceNodeIDs = elemNodes(localFaceNodeIDs);
        if hasVisited(elemFaces(i))
            tempFaceIDs = totalFaceIDs(:, elemFaces(i));
            
            newNodeCoors(:, i) = 0.5 * sum(elemNodeCoors(:, localFaceNodeIDs), 2);        
            if (faceTypes(i) == 0)
                newNodeIDs(i) = msh.faceNodes(2, tempFaceIDs(1));
            elseif (faceTypes(i) == 1)
                if (faceNodes(1, i) == faceNodeIDs(1) && faceNodes(2, i) == faceNodeIDs(2))
                    newNodeIDs(i) = msh.faceNodes(2, tempFaceIDs(1));
                elseif (faceNodes(1, i) == faceNodeIDs(2) && faceNodes(2, i) == faceNodeIDs(1))
                    newNodeIDs(i) = msh.faceNodes(2, tempFaceIDs(1));
                else
                    newNodeIDs(i) = msh.nNodes + 1;
                    msh.nodeCoor(:, newNodeIDs(i)) = newNodeCoors(:, i);
                    msh.nNodes = msh.nNodes + 1;    
                end
            else
                error('Wrong face type')
            end
            
            if (faceElems(1, i) == ID) 
                newFaceIDs(:, i) = tempFaceIDs;                
            else
                newFaceIDs(:, i) = flip(tempFaceIDs);                 
            end
        else
            newNodeIDs(i) = msh.nNodes + 1;           
            msh.nNodes    = msh.nNodes + 1;   
            if (faceTypes(i) == 0)
                newNodeCoors(:, i) = 0.5 * sum(elemNodeCoors(:, localFaceNodeIDs), 2);
                nodeID = newNodeIDs(i);
            elseif (faceTypes(i) == 1)
                newNodeCoors(:, i) = 0.5 * sum(elemNodeCoors(:, localFaceNodeIDs), 2);
                if (faceNodes(1, i) == faceNodeIDs(1) && faceNodes(2, i) == faceNodeIDs(2))
                    nodeID = newNodeIDs(i);
                elseif (faceNodes(1, i) == faceNodeIDs(2) && faceNodes(2, i) == faceNodeIDs(1))
                    nodeID = newNodeIDs(i);
                else
                    nodeID = msh.nNodes + 1;
                    msh.nodeCoor(:, nodeID) = 0.5 * sum(msh.nodeCoor(:, faceNodes(:, i)), 2);
                    msh.nNodes = msh.nNodes + 1; 
                end
            else
                newNodeCoors(:, i) = getBndNode(elemNodeCoors(:, localFaceNodeIDs(1)), elemNodeCoors(:, localFaceNodeIDs(2)), faceTypes(i));
                nodeID = newNodeIDs(i); 
            end
            msh.nodeCoor(:, newNodeIDs(i)) = newNodeCoors(:, i);
                    
            tempFaceIDs = msh.nFaces + [1; 2];
            totalFaceIDs(:, elemFaces(i)) = tempFaceIDs;
            msh.faceNodes(:, tempFaceIDs) = [faceNodes(1, i), nodeID; nodeID, faceNodes(2, i)];             
            if (faceTypes(i) == 0 || faceTypes(i) == 1)    
                if (faceElems(1, i) == ID)
                    newFaceIDs(:, i) = tempFaceIDs;
                    msh.faceElems(1, tempFaceIDs) = msh.nElems + getLocalSubelementIDs(i, true);
                    msh.faceElems(2, tempFaceIDs) = 4 * (faceElems(2, i) - 1) + getLocalSubelementIDs(faceNums(2, i), false);
                    msh.faceNums( 1, tempFaceIDs) = i;
                    msh.faceNums( 2, tempFaceIDs) = faceNums(2, i);
                else                   
                    newFaceIDs(:, i) = flip(tempFaceIDs);
                    msh.faceElems(1, tempFaceIDs) = 4 * (faceElems(1, i) - 1) + getLocalSubelementIDs(faceNums(1, i), true);
                    msh.faceElems(2, tempFaceIDs) = msh.nElems + getLocalSubelementIDs(i, false);
                    msh.faceNums( 1, tempFaceIDs) = faceNums(1, i);
                    msh.faceNums( 2, tempFaceIDs) = i;
                end                
            else 
                newFaceIDs(:, i) = tempFaceIDs;
                msh.faceElems(1, tempFaceIDs) = msh.nElems + getLocalSubelementIDs(i, true);
                msh.faceNums( 1, tempFaceIDs) = i;
            end
            if faceIsCurved(i)
                tempNormals = [newNodeCoors(2, i) - elemNodeCoors(2, localFaceNodeIDs(1)), elemNodeCoors(2, localFaceNodeIDs(2)) - newNodeCoors(2, i); elemNodeCoors(1, localFaceNodeIDs(1)) - newNodeCoors(1, i), newNodeCoors(1, i) - elemNodeCoors(1, localFaceNodeIDs(2))];
                tempSizes   = vecnorm(tempNormals);
                tempNormals = tempNormals ./ tempSizes;

                msh.faceNormalx(tempFaceIDs)  = tempNormals(1, :);
                msh.faceNormaly(tempFaceIDs)  = tempNormals(2, :);
                msh.faceSize(tempFaceIDs)     = tempSizes;
                msh.faceJac(tempFaceIDs)      = tempSizes / 2; 
                msh.faceIsCurved(tempFaceIDs) = true;
                msh.elemIsCurved(msh.faceElems(1, tempFaceIDs)) = true;                
            else
                msh.faceNormalx(tempFaceIDs) = faceNormalxs(i);
                msh.faceNormaly(tempFaceIDs) = faceNormalys(i);
                msh.faceSize(tempFaceIDs)    = faceSizes(i) / 2;
                msh.faceJac(tempFaceIDs)     = faceSizes(i) / 4;                 
            end
            msh.faceType(tempFaceIDs) = faceTypes(i);
            msh.nFaces = msh.nFaces + 2;    
            hasVisited(elemFaces(i)) = true;                 
        end
    end
    % interior nodes and four interior faces
    newNodeCoors(:, 5) = 0.25 * sum(newNodeCoors(:, 1 : 4), 2);
    newNodeIDs(5) = msh.nNodes + 1;
    msh.nodeCoor(:, newNodeIDs(5)) = newNodeCoors(:, 5);
    msh.nNodes = msh.nNodes + 1;
      
    newFaceIDs(:, 5 : 6) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [2, 2]);
    tempFaceIDs = msh.nFaces + [1; 2];
    msh.faceNodes(:, tempFaceIDs) = [newNodeIDs(3), newNodeIDs(5); newNodeIDs(5), newNodeIDs(4)];  
    msh.faceElems(1, tempFaceIDs) = msh.nElems + [1, 3];
    msh.faceElems(2, tempFaceIDs) = msh.nElems + [2, 4];
    msh.faceNums( 1, tempFaceIDs) = 2;
    msh.faceNums( 2, tempFaceIDs) = 1; 
    tempNormals = [newNodeCoors(2, 5) - newNodeCoors(2, 3), newNodeCoors(2, 4) - newNodeCoors(2, 5); newNodeCoors(1, 3) - newNodeCoors(1, 5), newNodeCoors(1, 5) - newNodeCoors(1, 4)];
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    msh.faceNormalx(tempFaceIDs) = tempNormals(1, :);
    msh.faceNormaly(tempFaceIDs) = tempNormals(2, :);
    msh.faceSize(tempFaceIDs)    = tempSizes;
    msh.faceJac(tempFaceIDs)     = tempSizes / 2;     
    tempFaceIDs = msh.nFaces + [3; 4];
    msh.faceNodes(:, tempFaceIDs) = [newNodeIDs(2), newNodeIDs(5); newNodeIDs(5), newNodeIDs(1)]; 
    msh.faceElems(1, tempFaceIDs) = msh.nElems + [2, 1];
    msh.faceElems(2, tempFaceIDs) = msh.nElems + [4, 3];
    msh.faceNums( 1, tempFaceIDs) = 4;
    msh.faceNums( 2, tempFaceIDs) = 3; 
    tempNormals = [newNodeCoors(2, 5) - newNodeCoors(2, 2), newNodeCoors(2, 1) - newNodeCoors(2, 5); newNodeCoors(1, 2) - newNodeCoors(1, 5), newNodeCoors(1, 5) - newNodeCoors(1, 1)];
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    msh.faceNormalx(tempFaceIDs) = tempNormals(1, :);
    msh.faceNormaly(tempFaceIDs) = tempNormals(2, :);
    msh.faceSize(tempFaceIDs)    = tempSizes;
    msh.faceJac(tempFaceIDs)     = tempSizes / 2;     
    msh.nFaces = msh.nFaces + 4;
    
    % element information
    tempElemIDs = msh.nElems + (1 : 4);
    elemNodeCoorx = [elemNodeCoors(1, 1),  newNodeCoors(1, 3),  newNodeCoors(1, 1),  newNodeCoors(1, 5);...
                      newNodeCoors(1, 3), elemNodeCoors(1, 2),  newNodeCoors(1, 5),  newNodeCoors(1, 2);...
                      newNodeCoors(1, 5),  newNodeCoors(1, 2),  newNodeCoors(1, 4), elemNodeCoors(1, 3);...
                      newNodeCoors(1, 1),  newNodeCoors(1, 5), elemNodeCoors(1, 4),  newNodeCoors(1, 4)]; 
                    
    elemNodeCoory = [elemNodeCoors(2, 1),  newNodeCoors(2, 3),  newNodeCoors(2, 1),  newNodeCoors(2, 5);...
                      newNodeCoors(2, 3), elemNodeCoors(2, 2),  newNodeCoors(2, 5),  newNodeCoors(2, 2);...
                      newNodeCoors(2, 5),  newNodeCoors(2, 2),  newNodeCoors(2, 4), elemNodeCoors(2, 3);...
                      newNodeCoors(2, 1),  newNodeCoors(2, 5), elemNodeCoors(2, 4),  newNodeCoors(2, 4)];
                     
    x    = quad.N     * elemNodeCoorx;
    xxi  = quad.N_xi  * elemNodeCoorx;
    xeta = quad.N_eta * elemNodeCoorx;
    y    = quad.N     * elemNodeCoory;
    yxi  = quad.N_xi  * elemNodeCoory;
    yeta = quad.N_eta * elemNodeCoory;     
    J    = xxi .* yeta - xeta .* yxi;
                
    si4 = 0.5 * ((elemNodeCoorx(1, :) - elemNodeCoorx(3, :)) .* (elemNodeCoory(2, :) - elemNodeCoory(4, :)) - (elemNodeCoorx(2, :) - elemNodeCoorx(4, :)) .* (elemNodeCoory(1, :) - elemNodeCoory(3, :)));  
    msh.elemSize(tempElemIDs)     = si4;                   
    msh.elemDiam(tempElemIDs)     = sqrt(si4);   
    msh.elemNodes(:, tempElemIDs) = [ elemNodes(1), newNodeIDs(3), newNodeIDs(1), newNodeIDs(5);...
                                     newNodeIDs(3),  elemNodes(2), newNodeIDs(5), newNodeIDs(2);...
                                     newNodeIDs(5), newNodeIDs(2), newNodeIDs(4),  elemNodes(3);...
                                     newNodeIDs(1), newNodeIDs(5),  elemNodes(4), newNodeIDs(4)];  
    msh.elemFaces(:, tempElemIDs) = [newFaceIDs(2, 1), newFaceIDs(1, 5), newFaceIDs(1, 1), newFaceIDs(2, 5);...
                                     newFaceIDs(1, 5), newFaceIDs(1, 2), newFaceIDs(2, 5), newFaceIDs(2, 2);...
                                     newFaceIDs(1, 3), newFaceIDs(2, 3), newFaceIDs(2, 6), newFaceIDs(1, 6);...
                                     newFaceIDs(2, 6), newFaceIDs(1, 6), newFaceIDs(2, 4), newFaceIDs(1, 4)];  
    msh.elemGPx(:, tempElemIDs)   = x;
    msh.elemGPy(:, tempElemIDs)   = y;
    msh.elemJac(:, tempElemIDs)   = J;
    msh.elemJxix(:, tempElemIDs)  = yeta; 
    msh.elemJxiy(:, tempElemIDs)  = -xeta; 
    msh.elemJetax(:, tempElemIDs) = -yxi; 
    msh.elemJetay(:, tempElemIDs) = xxi;                                    
    msh.nElems = msh.nElems + 4;    
    
    % inverse of mass matrix
    msh.massMatInv{tempElemIDs(1)} = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 1)));
    msh.massMatInv{tempElemIDs(2)} = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 2)));
    msh.massMatInv{tempElemIDs(3)} = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 3)));
    msh.massMatInv{tempElemIDs(4)} = inv(bs.phitw{1} * (bs.phi{1} .* J(:, 4)));    
end
msh.nodeCoor(:, msh.nNodes + 1 : end) = [];

% Update the IDs of internal and various boundary faces
isInternal = msh.faceType(msh.LFaces) == 0;
msh.intLFaces = msh.LFaces(isInternal);
msh.nIntLFaces = length(msh.intLFaces);
bndLFaces = msh.LFaces(~isInternal);
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = bndLFaces(msh.faceType(bndLFaces) == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (msh.type == 212)
    % Augment information of elements for adaptive mesh
    msh.elemLevel    = zeros(1, msh.nElems);
    msh.elemLID      = 1 : msh.nElems;
    msh.elemParent   = zeros(1, msh.nElems);
    msh.elemChildren = zeros(4, msh.nElems);
    
    % Augment information of faces for adaptive mesh
    msh.faceLevel    = zeros(1, msh.nFaces);
    msh.faceLID      = 1 : msh.nFaces;
    msh.faceParent   = zeros(1, msh.nFaces);
    msh.faceChildren = zeros(2, msh.nFaces);
end

function nodeIDs = getLocalFaceNodeIDs(faceID)

switch faceID
    case 1
        nodeIDs = [4, 1];
    case 2
        nodeIDs = [2, 3];
    case 3
        nodeIDs = [1, 2];
    case 4
        nodeIDs = [3, 4];
    otherwise
        error('Wrong face ID')
end
        
end

function subeIDs = getLocalSubelementIDs(faceID, isLeft)
    
switch faceID
    case 1
        subeIDs = [3, 1];
    case 2
        subeIDs = [2, 4];
    case 3
        subeIDs = [1, 2];
    case 4
        subeIDs = [4, 3];
    otherwise
        error('Wrong face ID')
end
if ~isLeft
    subeIDs = flip(subeIDs);
end

end

end


