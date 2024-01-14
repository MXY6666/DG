% Copyright (C) Xiangyi Meng.

function msh = refineHexaMeshGlobally(msh0, getFaceNode, getEdgeNode, quad, bs)

if (msh0.type ~= 311) && (msh0.type ~= 312)
    error('Wrong mesh type')
end
if (msh0.nElems ~= msh0.nLElems)
    error('All the elements of the background mesh should be leaf elements')
end

if (nargin < 2) || isempty(getFaceNode)
    getFaceNode = @(p1, p2, p3, p4, faceType)0.25 * (p1 + p2 + p3 + p4);
end

if (nargin < 3) || isempty(getEdgeNode)
    getEdgeNode = @(p1, p2, edgeType)0.5 * (p1 + p2);
end

if (nargin < 4) || isempty(quad)
    quad = GaussQuadratureRule_hexa([2, 2, 2], 301);
end

if (nargin < 5) || isempty(bs)
    quad2 = GaussQuadratureRule_square([2, 2], 201);
    bs    = setBasisFunctionSet_cube(quad2, quad, 301);
end
%**************************************************************************
n12 = msh0.N(1) * msh0.N(2); n23 = msh0.N(2) * msh0.N(3); n31 = msh0.N(3) * msh0.N(1);
msh          = getEmptyMesh;
msh.dm       = msh0.dm;
msh.N        = 2 * msh0.N;
msh.type     = msh0.type;
msh.maxLevel = msh0.maxLevel;
msh.bndTypes = msh0.bndTypes;
msh.nElems   = 0;
msh.nFaces   = 0;
if (msh0.type == 312)
    msh.nEdges = 0;
end
msh.nNodes   = msh0.nNodes;
msh.nLElems  = 8 * msh0.nElems; 
msh.nLFaces  = 12 * msh0.nElems + 4 * msh0.nFaces; 
msh.LElems   = 1 : msh.nLElems;
msh.LFaces   = 1 : msh.nLFaces;
msh.nGPs     = msh0.nGPs;
%**************************************************************************
si = 8 * msh0.nElems;
msh.elemSize     = zeros(1, si);
msh.elemDiam     = zeros(1, si);
msh.elemNodes    = zeros(8, si);
msh.elemEdges    = zeros(12, si);
msh.elemFaces    = zeros(6, si);
msh.elemGPx      = zeros(quad.np, si);
msh.elemGPy      = zeros(quad.np, si);
msh.elemGPz      = zeros(quad.np, si);
msh.elemJac      = zeros(quad.np, si);
msh.elemJxix     = zeros(quad.np, si); 
msh.elemJxiy     = zeros(quad.np, si); 
msh.elemJxiz     = zeros(quad.np, si);
msh.elemJetax    = zeros(quad.np, si); 
msh.elemJetay    = zeros(quad.np, si); 
msh.elemJetaz    = zeros(quad.np, si); 
msh.elemJzetax   = zeros(quad.np, si); 
msh.elemJzetay   = zeros(quad.np, si); 
msh.elemJzetaz   = zeros(quad.np, si); 
msh.elemIsCurved = false(1, si);
si = 12 * msh0.nElems + 4 * msh0.nFaces;
msh.faceNormalx  = zeros(1, si); 
msh.faceNormaly  = zeros(1, si);
msh.faceNormalz  = zeros(1, si);
msh.faceSize     = zeros(1, si);
msh.faceType     = zeros(1, si);
msh.faceNodes    = zeros(4, si);
msh.faceEdges    = zeros(4, si);
msh.faceElems    = zeros(2, si); 
msh.faceNums     = zeros(2, si); 
msh.faceR2L      = zeros(1, si);
msh.faceJac      = zeros(1, si); 
msh.faceIsCurved = false(1, si);
if (msh0.type == 311)
    v1 = [1, 4, 5, 8, 1, 2, 5, 6, 1, 2, 4, 3];
    v2 = [2, 3, 6, 7, 4, 3, 8, 7, 5, 6, 8, 7];
    Q = sparse(msh0.elemNodes(v1, :), msh0.elemNodes(v2, :), 1, msh0.nNodes, msh0.nNodes);
    [I, J] = find(Q | Q');
    K = find(I >= J);
    I = I(K);
    J = J(K);
    ned0 = length(I);
    Q = sparse(I, J, 1 : ned0, msh0.nNodes, msh0.nNodes);
    
    elemEdgeIDs = sub2ind([msh0.nNodes, msh0.nNodes], max(msh0.elemNodes(v1, :), msh0.elemNodes(v2, :)), min(msh0.elemNodes(v1, :), msh0.elemNodes(v2, :)));
    elemEdgeIDs = reshape(nonzeros(Q(elemEdgeIDs)), [12, msh0.nElems]);
    faceEdgeIDs = sub2ind([msh0.nNodes, msh0.nNodes], max(msh0.faceNodes([1, 2, 1, 4], :), msh0.faceNodes([4, 3, 2, 3], :)), min(msh0.faceNodes([1, 2, 1, 4], :), msh0.faceNodes([4, 3, 2, 3], :)));
    faceEdgeIDs = reshape(nonzeros(Q(faceEdgeIDs)), [4, msh0.nFaces]);  
else
    si = 6 * msh0.nElems + 4 * msh0.nFaces + 2 * msh0.nEdges + 4 * (n12 + n23 + n31);
    msh.edgeType     = zeros(1, si); 
    msh.edgeNodes    = zeros(2, si); 
    msh.edgeNElems   = zeros(1, si);  
    ned0 = msh0.nEdges;
end
si = msh0.nNodes + msh0.nElems + msh0.nFaces + ned0 + n12 + n23 + n31;
msh.nodeCoor = zeros(3, si);
msh.nodeCoor(:, 1 : msh0.nNodes) = msh0.nodeCoor;
%**************************************************************************
faceHasVisited = false(1, msh0.nFaces);
edgeHasVisited = false(1, ned0);
newFaceIDs     = zeros(4, msh0.nFaces);
newFaceNodeIDs = zeros(1, msh0.nFaces);
newEdgeNodeIDs = zeros(1, ned0);
ind = [1,2,4,5,10,11,13,14;2,3,5,6,11,12,14,15;5,6,8,9,14,15,17,18;4,5,7,8,13,14,16,17;10,11,13,14,19,20,22,23;11,12,14,15,20,21,23,24;14,15,17,18,23,24,26,27;13,14,16,17,22,23,25,26];
faceIDs1  = zeros(3, 2, 2);
faceIDs2  = zeros(2, 3, 2);
faceIDs3  = zeros(2, 2, 3);
nodeIDs   = zeros(3, 3, 3);
nodeCoors = zeros(3, 27);
if (msh0.type == 312)
    newFaceEdgeIDs = zeros(4, msh0.nFaces);
    newEdgeIDs     = zeros(2, ned0);    
    edgeIDs1       = zeros(2, 3, 3);
    edgeIDs2       = zeros(3, 2, 3);
    edgeIDs3       = zeros(3, 3, 2);     
end

if (msh0.type == 311)
    for ID = 1 : msh0.nElems
        % current loop element information
        elemNodes     = msh0.elemNodes(:, ID);
        elemEdges     = elemEdgeIDs(:, ID);
        elemFaces     = msh0.elemFaces(:, ID);
        elemNodeCoors = msh0.nodeCoor(:, elemNodes);
        
        nodeIDs([1, 3, 9, 7, 19, 21, 27, 25])      = elemNodes;
        nodeCoors(:, [1, 3, 9, 7, 19, 21, 27, 25]) = elemNodeCoors;
      
        % cration of elements
        elemIDs = reshape(msh.nElems + 1 : msh.nElems + 8, [2, 2, 2]);

        % Deal with the edges
        [nodeIDs(2),  nodeCoors(:, 2)]  = refineEdge1(elemEdges(1));  % 1st  edge
        [nodeIDs(8),  nodeCoors(:, 8)]  = refineEdge1(elemEdges(2));  % 2nd  edge
        [nodeIDs(20), nodeCoors(:, 20)] = refineEdge1(elemEdges(3));  % 3rd  edge
        [nodeIDs(26), nodeCoors(:, 26)] = refineEdge1(elemEdges(4));  % 4th  edge
        [nodeIDs(4),  nodeCoors(:, 4)]  = refineEdge1(elemEdges(5));  % 5th  edge
        [nodeIDs(6),  nodeCoors(:, 6)]  = refineEdge1(elemEdges(6));  % 6th  edge
        [nodeIDs(22), nodeCoors(:, 22)] = refineEdge1(elemEdges(7));  % 7th  edge
        [nodeIDs(24), nodeCoors(:, 24)] = refineEdge1(elemEdges(8));  % 8th  edge
        [nodeIDs(10), nodeCoors(:, 10)] = refineEdge1(elemEdges(9));  % 9th  edge
        [nodeIDs(12), nodeCoors(:, 12)] = refineEdge1(elemEdges(10)); % 10th edge
        [nodeIDs(16), nodeCoors(:, 16)] = refineEdge1(elemEdges(11)); % 11th edge
        [nodeIDs(18), nodeCoors(:, 18)] = refineEdge1(elemEdges(12)); % 12th edge
        
        % Deal with the faces
        % 1st face
        [nodeIDs(13), nodeCoors(:, 13), faceIDs] = refineFace1(ID, elemFaces(1), elemNodes([1, 4, 8, 5]), elemNodeCoors(:, [1, 4, 8, 5]));
        faceIDs1(1, :, :) = reshape(faceIDs, [1, 2, 2]);
        % 2nd face
        [nodeIDs(15), nodeCoors(:, 15), faceIDs] = refineFace1(ID, elemFaces(2), elemNodes([2, 3, 7, 6]), elemNodeCoors(:, [2, 3, 7, 6]));
        faceIDs1(3, :, :) = reshape(faceIDs, [1, 2, 2]);    
        % 3rd face
        [nodeIDs(11), nodeCoors(:, 11), faceIDs] = refineFace1(ID, elemFaces(3), elemNodes([1, 5, 6, 2]), elemNodeCoors(:, [1, 5, 6, 2]));
        faceIDs2(:, 1, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);       
        % 4th face
        [nodeIDs(17), nodeCoors(:, 17), faceIDs] = refineFace1(ID, elemFaces(4), elemNodes([4, 8, 7, 3]), elemNodeCoors(:, [4, 8, 7, 3]));
        faceIDs2(:, 3, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);      
        % 5th face
        [nodeIDs(5), nodeCoors(:, 5), faceIDs] = refineFace1(ID, elemFaces(5), elemNodes(1 : 4), elemNodeCoors(:, 1 : 4));
        faceIDs3(:, :, 1) = reshape(faceIDs, [2, 2]);      
        % 6th face
        [nodeIDs(23), nodeCoors(:, 23), faceIDs] = refineFace1(ID, elemFaces(6), elemNodes(5 : 8), elemNodeCoors(:, 5 : 8));
        faceIDs3(:, :, 3) = reshape(faceIDs, [2, 2]); 
        
        % interior nodes, edges and faces
        nodeIDs(14) = msh.nNodes + 1;
        nodeCoors(:, 14) = sum(nodeCoors(:, [5, 11, 13, 15, 17, 23]), 2) / 6;
        msh.nodeCoor(:, nodeIDs(14)) = nodeCoors(:, 14);
        msh.nNodes = msh.nNodes + 1;
   
        faceIDs1(2, :, :) = reshape(msh.nFaces + 1 : msh.nFaces + 4, [1, 2, 2]);
        msh.faceNodes(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 1 : 2, 1 : 2), [1, 4]);
        msh.faceNodes(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 2 : 3, 1 : 2), [1, 4]);
        msh.faceNodes(3, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 2 : 3, 2 : 3), [1, 4]);
        msh.faceNodes(4, msh.nFaces + 1 : msh.nFaces + 4) = reshape(nodeIDs(2, 1 : 2, 2 : 3), [1, 4]);
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(1, :, :), [1, 4]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(2, :, :), [1, 4]);
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
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = elemIDs([1, 5, 2, 6]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = elemIDs([3, 7, 4, 8]);
        msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 4;
        msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 3;    
        msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4)      = 1;
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
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(:, :, 1), [1, 4]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(:, :, 2), [1, 4]);
        msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 6;
        msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 5;    
        msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4)      = 1;
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
        elemNodeCoorx = reshape(nodeCoors(1, ind), [8, 8]);
        elemNodeCoory = reshape(nodeCoors(2, ind), [8, 8]);
        elemNodeCoorz = reshape(nodeCoors(3, ind), [8, 8]);

        x     = quad.N      * elemNodeCoorx;
        xxi   = quad.N_xi   * elemNodeCoorx;
        xeta  = quad.N_eta  * elemNodeCoorx;
        xzeta = quad.N_zeta * elemNodeCoorx;
        y     = quad.N      * elemNodeCoory;
        yxi   = quad.N_xi   * elemNodeCoory;
        yeta  = quad.N_eta  * elemNodeCoory;   
        yzeta = quad.N_zeta * elemNodeCoory;  
        z     = quad.N      * elemNodeCoorz;
        zxi   = quad.N_xi   * elemNodeCoorz;
        zeta  = quad.N_eta  * elemNodeCoorz;
        zzeta = quad.N_zeta * elemNodeCoorz;  
        Jac   = xxi .* (yeta .* zzeta - yzeta .* zeta) + xeta .* (yzeta .* zxi - yxi .* zzeta) + xzeta .* (yxi .* zeta - yeta .* zxi);

        si = 1 / 6 * dot(nodeCoors(:, ind(7, :)) - nodeCoors(:, ind(1, :)), mycross3(nodeCoors(:, ind(2, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(3, :)) - nodeCoors(:, ind(6, :))) + mycross3(nodeCoors(:, ind(4, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(8, :)) - nodeCoors(:, ind(3, :))) + mycross3(nodeCoors(:, ind(5, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(6, :)) - nodeCoors(:, ind(8, :))));
        msh.elemSize(elemIDs)      =  si;
        msh.elemDiam(elemIDs)      =  si.^(1/3);
        msh.elemNodes(1, elemIDs)  = reshape(nodeIDs(1 : 2, 1 : 2, 1 : 2), [1, 8]);
        msh.elemNodes(2, elemIDs)  = reshape(nodeIDs(2 : 3, 1 : 2, 1 : 2), [1, 8]);
        msh.elemNodes(3, elemIDs)  = reshape(nodeIDs(2 : 3, 2 : 3, 1 : 2), [1, 8]);
        msh.elemNodes(4, elemIDs)  = reshape(nodeIDs(1 : 2, 2 : 3, 1 : 2), [1, 8]);
        msh.elemNodes(5, elemIDs)  = reshape(nodeIDs(1 : 2, 1 : 2, 2 : 3), [1, 8]);
        msh.elemNodes(6, elemIDs)  = reshape(nodeIDs(2 : 3, 1 : 2, 2 : 3), [1, 8]);
        msh.elemNodes(7, elemIDs)  = reshape(nodeIDs(2 : 3, 2 : 3, 2 : 3), [1, 8]);
        msh.elemNodes(8, elemIDs)  = reshape(nodeIDs(1 : 2, 2 : 3, 2 : 3), [1, 8]);   
        msh.elemFaces(1, elemIDs)  = reshape(faceIDs1(1 : 2, :, :), [1, 8]);
        msh.elemFaces(2, elemIDs)  = reshape(faceIDs1(2 : 3, :, :), [1, 8]);
        msh.elemFaces(3, elemIDs)  = reshape(faceIDs2(:, 1 : 2, :), [1, 8]);
        msh.elemFaces(4, elemIDs)  = reshape(faceIDs2(:, 2 : 3, :), [1, 8]);
        msh.elemFaces(5, elemIDs)  = reshape(faceIDs3(:, :, 1 : 2), [1, 8]);
        msh.elemFaces(6, elemIDs)  = reshape(faceIDs3(:, :, 2 : 3), [1, 8]); 
        msh.elemGPx(:, elemIDs)    = x;
        msh.elemGPy(:, elemIDs)    = y;
        msh.elemGPz(:, elemIDs)    = z;
        msh.elemJac(:, elemIDs)    = Jac;
        msh.elemJxix(:, elemIDs)   = yeta .* zzeta - yzeta .* zeta;
        msh.elemJxiy(:, elemIDs)   = zeta .* xzeta - zzeta .* xeta;
        msh.elemJxiz(:, elemIDs)   = xeta .* yzeta - xzeta .* yeta;
        msh.elemJetax(:, elemIDs)  = yzeta .* zxi - yxi .* zzeta;
        msh.elemJetay(:, elemIDs)  = zzeta .* xxi - zxi .* xzeta;
        msh.elemJetaz(:, elemIDs)  = xzeta .* yxi - xxi .* yzeta;
        msh.elemJzetax(:, elemIDs) = yxi .* zeta - yeta .* zxi;
        msh.elemJzetay(:, elemIDs) = zxi .* xeta - zeta .* xxi;
        msh.elemJzetaz(:, elemIDs) = xxi .* yeta - xeta .* yxi;
        msh.nElems = msh.nElems + 8;   
        
        % inverse of mass matrix
        msh.massMatInv{elemIDs(1)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 1)));
        msh.massMatInv{elemIDs(2)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 2)));
        msh.massMatInv{elemIDs(3)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 3)));
        msh.massMatInv{elemIDs(4)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 4)));
        msh.massMatInv{elemIDs(5)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 5)));
        msh.massMatInv{elemIDs(6)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 6)));
        msh.massMatInv{elemIDs(7)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 7)));
        msh.massMatInv{elemIDs(8)} = inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 8)));        
    end
else
    for ID = 1 : msh0.nElems
        % current loop element information
        elemNodes     = msh0.elemNodes(:, ID);
        elemEdges     = msh0.elemEdges(:, ID);
        elemFaces     = msh0.elemFaces(:, ID);
        elemNodeCoors = msh0.nodeCoor(:, elemNodes);
        elemEdgeNodes = msh0.edgeNodes(:, elemEdges);
        
        nodeIDs([1, 3, 9, 7, 19, 21, 27, 25])      = elemNodes;
        nodeCoors(:, [1, 3, 9, 7, 19, 21, 27, 25]) = elemNodeCoors;
        
        % cration of elements
        elemIDs = reshape(msh.nElems + 1 : msh.nElems + 8, [2, 2, 2]);

        % Deal with the edges
        % 1st edge
        [nodeIDs(2), nodeCoors(:, 2), edgeIDs] = refineEdge2(elemEdges(1));
        if (elemEdgeNodes(1, 1) == nodeIDs(1))
            edgeIDs1(:, 1, 1) = edgeIDs;
        else
            edgeIDs1(:, 1, 1) = flip(edgeIDs);        
        end
        % 2nd edge
        [nodeIDs(8), nodeCoors(:, 8), edgeIDs] = refineEdge2(elemEdges(2));
        if (elemEdgeNodes(1, 2) == nodeIDs(7))
            edgeIDs1(:, 3, 1) = edgeIDs;
        else
            edgeIDs1(:, 3, 1) = flip(edgeIDs);        
        end       
        % 3rd edge
        [nodeIDs(20), nodeCoors(:, 20), edgeIDs] = refineEdge2(elemEdges(3));
        if (elemEdgeNodes(1, 3) == nodeIDs(19))
            edgeIDs1(:, 1, 3) = edgeIDs;
        else
            edgeIDs1(:, 1, 3) = flip(edgeIDs);        
        end      
        % 4th edge
        [nodeIDs(26), nodeCoors(:, 26), edgeIDs] = refineEdge2(elemEdges(4));
        if (elemEdgeNodes(1, 4) == nodeIDs(25))
            edgeIDs1(:, 3, 3) = edgeIDs;
        else
            edgeIDs1(:, 3, 3) = flip(edgeIDs);        
        end       
        % 5th edge
        [nodeIDs(4), nodeCoors(:, 4), edgeIDs] = refineEdge2(elemEdges(5));
        if (elemEdgeNodes(1, 5) == nodeIDs(1))
            edgeIDs2(1, :, 1) = edgeIDs;
        else
            edgeIDs2(1, :, 1) = flip(edgeIDs);        
        end     
        % 6th edge
        [nodeIDs(6), nodeCoors(:, 6), edgeIDs] = refineEdge2(elemEdges(6));
        if (elemEdgeNodes(1, 6) == nodeIDs(3))
            edgeIDs2(3, :, 1) = edgeIDs;
        else
            edgeIDs2(3, :, 1) = flip(edgeIDs);        
        end   
        % 7th edge
        [nodeIDs(22), nodeCoors(:, 22), edgeIDs] = refineEdge2(elemEdges(7));
        if (elemEdgeNodes(1, 7) == nodeIDs(19))
            edgeIDs2(1, :, 3) = edgeIDs;
        else
            edgeIDs2(1, :, 3) = flip(edgeIDs);        
        end  
        % 8th edge
        [nodeIDs(24), nodeCoors(:, 24), edgeIDs] = refineEdge2(elemEdges(8));
        if (elemEdgeNodes(1, 8) == nodeIDs(21))
            edgeIDs2(3, :, 3) = edgeIDs;
        else
            edgeIDs2(3, :, 3) = flip(edgeIDs);        
        end    
        % 9th edge
        [nodeIDs(10), nodeCoors(:, 10), edgeIDs] = refineEdge2(elemEdges(9));
        if (elemEdgeNodes(1, 9) == nodeIDs(1))
            edgeIDs3(1, 1, :) = edgeIDs;
        else
            edgeIDs3(1, 1, :) = flip(edgeIDs);        
        end 
        % 10th edge
        [nodeIDs(12), nodeCoors(:, 12), edgeIDs] = refineEdge2(elemEdges(10));
        if (elemEdgeNodes(1, 10) == nodeIDs(3))
            edgeIDs3(3, 1, :) = edgeIDs;
        else
            edgeIDs3(3, 1, :) = flip(edgeIDs);        
        end 
        % 11th edge
        [nodeIDs(16), nodeCoors(:, 16), edgeIDs] = refineEdge2(elemEdges(11));
        if (elemEdgeNodes(1, 11) == nodeIDs(7))
            edgeIDs3(1, 3, :) = edgeIDs;
        else
            edgeIDs3(1, 3, :) = flip(edgeIDs);        
        end    
        % 12th edge
        [nodeIDs(18), nodeCoors(:, 18), edgeIDs] = refineEdge2(elemEdges(12));
        if (elemEdgeNodes(1, 12) == nodeIDs(9))
            edgeIDs3(3, 3, :) = edgeIDs;
        else
            edgeIDs3(3, 3, :) = flip(edgeIDs);        
        end

        % Deal with the faces
        % 1st face
        elemFaceNodes = reshape(nodeIDs(1, :, :), [9, 1]);
        [nodeIDs(13), nodeCoors(:, 13), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(1), elemFaceNodes, elemNodeCoors(:, [1, 4, 8, 5]));
        edgeIDs3(1, 2, :) = edgeIDs(1 : 2);
        edgeIDs2(1, :, 2) = edgeIDs(3 : 4);
        faceIDs1(1, :, :) = reshape(faceIDs, [1, 2, 2]);   
        % 2nd face
        elemFaceNodes = reshape(nodeIDs(3, :, :), [9, 1]);
        [nodeIDs(15), nodeCoors(:, 15), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(2), elemFaceNodes, elemNodeCoors(:, [2, 3, 7, 6]));
        edgeIDs3(3, 2, :) = edgeIDs(1 : 2);
        edgeIDs2(3, :, 2) = edgeIDs(3 : 4);
        faceIDs1(3, :, :) = reshape(faceIDs, [1, 2, 2]);    
        % 3rd face
        elemFaceNodes = nodeIDs([1, 10, 19, 2, 11, 20, 3, 12, 21])';
        [nodeIDs(11), nodeCoors(:, 11), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(3), elemFaceNodes, elemNodeCoors(:, [1, 5, 6, 2]));
        edgeIDs1(:, 1, 2) = edgeIDs(1 : 2);
        edgeIDs3(2, 1, :) = edgeIDs(3 : 4);
        faceIDs2(:, 1, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);       
        % 4th face
        elemFaceNodes = nodeIDs([7, 16, 25, 8, 17, 26, 9, 18, 27])';
        [nodeIDs(17), nodeCoors(:, 17), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(4), elemFaceNodes, elemNodeCoors(:, [4, 8, 7, 3]));
        edgeIDs1(:, 3, 2) = edgeIDs(1 : 2);
        edgeIDs3(2, 3, :) = edgeIDs(3 : 4);
        faceIDs2(:, 3, :) = reshape(faceIDs([1, 3, 2, 4]), [2, 1, 2]);      
        % 5th face
        elemFaceNodes = reshape(nodeIDs(:, :, 1), [9, 1]);
        [nodeIDs(5), nodeCoors(:, 5), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(5), elemFaceNodes, elemNodeCoors(:, 1 : 4));
        edgeIDs2(2, :, 1) = edgeIDs(1 : 2);
        edgeIDs1(:, 2, 1) = edgeIDs(3 : 4);
        faceIDs3(:, :, 1) = reshape(faceIDs, [2, 2]);      
        % 6th face
        elemFaceNodes = reshape(nodeIDs(:, :, 3), [9, 1]);
        [nodeIDs(23), nodeCoors(:, 23), edgeIDs, faceIDs] = refineFace2(ID, elemFaces(6), elemFaceNodes, elemNodeCoors(:, 5 : 8));
        edgeIDs2(2, :, 3) = edgeIDs(1 : 2);
        edgeIDs1(:, 2, 3) = edgeIDs(3 : 4);
        faceIDs3(:, :, 3) = reshape(faceIDs, [2, 2]); 

        % interior nodes, edges and faces
        nodeIDs(14) = msh.nNodes + 1;
        nodeCoors(:, 14) = sum(nodeCoors(:, [5, 11, 13, 15, 17, 23]), 2) / 6;
        msh.nodeCoor(:, nodeIDs(14)) = nodeCoors(:, 14);
        msh.nNodes = msh.nNodes + 1;

        edgeIDs1(:, 2, 2) = msh.nEdges + 1 : msh.nEdges + 2;
        msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2) = [nodeIDs(13), nodeIDs(14); nodeIDs(14), nodeIDs(15)];
        msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)   = 4;  
        msh.nEdges = msh.nEdges + 2;   

        edgeIDs2(2, :, 2) = msh.nEdges + 1 : msh.nEdges + 2;
        msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2) = [nodeIDs(11), nodeIDs(14); nodeIDs(14), nodeIDs(17)];
        msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)   = 4;  
        msh.nEdges = msh.nEdges + 2;   

        edgeIDs3(2, 2, :) = msh.nEdges + 1 : msh.nEdges + 2;
        msh.edgeNodes(:, msh.nEdges + 1 : msh.nEdges + 2) = [nodeIDs(5), nodeIDs(14); nodeIDs(14), nodeIDs(23)];
        msh.edgeNElems(msh.nEdges + 1 : msh.nEdges + 2)   = 4;  
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
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(1, :, :), [1, 4]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(2, :, :), [1, 4]);
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
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = elemIDs([1, 5, 2, 6]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = elemIDs([3, 7, 4, 8]);
        msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 4;
        msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 3;    
        msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4)      = 1;
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
        msh.faceElems(1, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(:, :, 1), [1, 4]);
        msh.faceElems(2, msh.nFaces + 1 : msh.nFaces + 4) = reshape(elemIDs(:, :, 2), [1, 4]);
        msh.faceNums( 1, msh.nFaces + 1 : msh.nFaces + 4) = 6;
        msh.faceNums( 2, msh.nFaces + 1 : msh.nFaces + 4) = 5;    
        msh.faceR2L(msh.nFaces + 1 : msh.nFaces + 4)      = 1;
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
        elemNodeCoorx = reshape(nodeCoors(1, ind), [8, 8]);
        elemNodeCoory = reshape(nodeCoors(2, ind), [8, 8]);
        elemNodeCoorz = reshape(nodeCoors(3, ind), [8, 8]);

        x     = quad.N      * elemNodeCoorx;
        xxi   = quad.N_xi   * elemNodeCoorx;
        xeta  = quad.N_eta  * elemNodeCoorx;
        xzeta = quad.N_zeta * elemNodeCoorx;
        y     = quad.N      * elemNodeCoory;
        yxi   = quad.N_xi   * elemNodeCoory;
        yeta  = quad.N_eta  * elemNodeCoory;   
        yzeta = quad.N_zeta * elemNodeCoory;  
        z     = quad.N      * elemNodeCoorz;
        zxi   = quad.N_xi   * elemNodeCoorz;
        zeta  = quad.N_eta  * elemNodeCoorz;
        zzeta = quad.N_zeta * elemNodeCoorz;  
        Jac   = xxi .* (yeta .* zzeta - yzeta .* zeta) + xeta .* (yzeta .* zxi - yxi .* zzeta) + xzeta .* (yxi .* zeta - yeta .* zxi);

        si = 1 / 6 * dot(nodeCoors(:, ind(7, :)) - nodeCoors(:, ind(1, :)), mycross3(nodeCoors(:, ind(2, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(3, :)) - nodeCoors(:, ind(6, :))) + mycross3(nodeCoors(:, ind(4, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(8, :)) - nodeCoors(:, ind(3, :))) + mycross3(nodeCoors(:, ind(5, :)) - nodeCoors(:, ind(1, :)), nodeCoors(:, ind(6, :)) - nodeCoors(:, ind(8, :))));
        msh.elemSize(elemIDs)      = si;
        msh.elemDiam(elemIDs)      = si.^(1/3);
        msh.elemNodes(1, elemIDs)  = reshape(nodeIDs(1 : 2, 1 : 2, 1 : 2), [1, 8]);
        msh.elemNodes(2, elemIDs)  = reshape(nodeIDs(2 : 3, 1 : 2, 1 : 2), [1, 8]);
        msh.elemNodes(3, elemIDs)  = reshape(nodeIDs(2 : 3, 2 : 3, 1 : 2), [1, 8]);
        msh.elemNodes(4, elemIDs)  = reshape(nodeIDs(1 : 2, 2 : 3, 1 : 2), [1, 8]);
        msh.elemNodes(5, elemIDs)  = reshape(nodeIDs(1 : 2, 1 : 2, 2 : 3), [1, 8]);
        msh.elemNodes(6, elemIDs)  = reshape(nodeIDs(2 : 3, 1 : 2, 2 : 3), [1, 8]);
        msh.elemNodes(7, elemIDs)  = reshape(nodeIDs(2 : 3, 2 : 3, 2 : 3), [1, 8]);
        msh.elemNodes(8, elemIDs)  = reshape(nodeIDs(1 : 2, 2 : 3, 2 : 3), [1, 8]);
        msh.elemEdges(1, elemIDs)  = reshape(edgeIDs1(:, 1 : 2, 1 : 2), [1, 8]);
        msh.elemEdges(2, elemIDs)  = reshape(edgeIDs1(:, 2 : 3, 1 : 2), [1, 8]);
        msh.elemEdges(3, elemIDs)  = reshape(edgeIDs1(:, 1 : 2, 2 : 3), [1, 8]);
        msh.elemEdges(4, elemIDs)  = reshape(edgeIDs1(:, 2 : 3, 2 : 3), [1, 8]);  
        msh.elemEdges(5, elemIDs)  = reshape(edgeIDs2(1 : 2, :, 1 : 2), [1, 8]);
        msh.elemEdges(6, elemIDs)  = reshape(edgeIDs2(2 : 3, :, 1 : 2), [1, 8]);
        msh.elemEdges(7, elemIDs)  = reshape(edgeIDs2(1 : 2, :, 2 : 3), [1, 8]);
        msh.elemEdges(8, elemIDs)  = reshape(edgeIDs2(2 : 3, :, 2 : 3), [1, 8]);    
        msh.elemEdges(9, elemIDs)  = reshape(edgeIDs3(1 : 2, 1 : 2, :), [1, 8]);
        msh.elemEdges(10, elemIDs) = reshape(edgeIDs3(2 : 3, 1 : 2, :), [1, 8]);
        msh.elemEdges(11, elemIDs) = reshape(edgeIDs3(1 : 2, 2 : 3, :), [1, 8]);
        msh.elemEdges(12, elemIDs) = reshape(edgeIDs3(2 : 3, 2 : 3, :), [1, 8]);         
        msh.elemFaces(1, elemIDs)  = reshape(faceIDs1(1 : 2, :, :), [1, 8]);
        msh.elemFaces(2, elemIDs)  = reshape(faceIDs1(2 : 3, :, :), [1, 8]);
        msh.elemFaces(3, elemIDs)  = reshape(faceIDs2(:, 1 : 2, :), [1, 8]);
        msh.elemFaces(4, elemIDs)  = reshape(faceIDs2(:, 2 : 3, :), [1, 8]);
        msh.elemFaces(5, elemIDs)  = reshape(faceIDs3(:, :, 1 : 2), [1, 8]);
        msh.elemFaces(6, elemIDs)  = reshape(faceIDs3(:, :, 2 : 3), [1, 8]); 
        msh.elemGPx(:, elemIDs)    = x;
        msh.elemGPy(:, elemIDs)    = y;
        msh.elemGPz(:, elemIDs)    = z;
        msh.elemJac(:, elemIDs)    = Jac;
        msh.elemJxix(:, elemIDs)   = yeta .* zzeta - yzeta .* zeta;
        msh.elemJxiy(:, elemIDs)   = zeta .* xzeta - zzeta .* xeta;
        msh.elemJxiz(:, elemIDs)   = xeta .* yzeta - xzeta .* yeta;
        msh.elemJetax(:, elemIDs)  = yzeta .* zxi - yxi .* zzeta;
        msh.elemJetay(:, elemIDs)  = zzeta .* xxi - zxi .* xzeta;
        msh.elemJetaz(:, elemIDs)  = xzeta .* yxi - xxi .* yzeta;
        msh.elemJzetax(:, elemIDs) = yxi .* zeta - yeta .* zxi;
        msh.elemJzetay(:, elemIDs) = zxi .* xeta - zeta .* xxi;
        msh.elemJzetaz(:, elemIDs) = xxi .* yeta - xeta .* yxi;
        msh.nElems = msh.nElems + 8;   
        
        % inverse of mass matrix
        msh.massMatInv{elemIDs(1)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 1))), 1.0e-10);
        msh.massMatInv{elemIDs(2)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 2))), 1.0e-10);
        msh.massMatInv{elemIDs(3)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 3))), 1.0e-10);
        msh.massMatInv{elemIDs(4)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 4))), 1.0e-10);
        msh.massMatInv{elemIDs(5)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 5))), 1.0e-10);
        msh.massMatInv{elemIDs(6)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 6))), 1.0e-10);
        msh.massMatInv{elemIDs(7)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 7))), 1.0e-10);
        msh.massMatInv{elemIDs(8)} = mychop(inv(bs.phitw{1} * (bs.phi{1} .* Jac(:, 8))), 1.0e-10);             
    end  
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
if (msh0.type == 312)
    % Augment information of elements for adaptive mesh
    msh.elemLevel    = zeros(1, msh.nElems);
    msh.elemLID      = 1 : msh.nElems;
    msh.elemParent   = zeros(1, msh.nElems);
    msh.elemChildren = zeros(8, msh.nElems);
    
    % Augment information of faces for adaptive mesh
    msh.faceLevel    = zeros(1, msh.nFaces);
    msh.faceLID      = 1 : msh.nFaces;
    msh.faceParent   = zeros(1, msh.nFaces);
    msh.faceChildren = zeros(4, msh.nFaces);
    
    % Augment information of edges for adaptive mesh
    msh.edgeParent = zeros(1, msh.nEdges);
    msh.edgeChildren = zeros(2, msh.nEdges);
    
    msh.edgeType(msh.nEdges + 1 : end)     = []; 
    msh.edgeNodes(:, msh.nEdges + 1 : end) = []; 
    msh.edgeNElems(msh.nEdges + 1 : end)   = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine an edge
function [nodeID, nodeCoor] = refineEdge1(edgeID)
    
if edgeHasVisited(edgeID)    
    nodeID   = newEdgeNodeIDs(edgeID);
    nodeCoor = msh.nodeCoor(:, nodeID);
else
    % Warning : This may not be true for edges who lie on a curved faces
    nodeID = msh.nNodes + 1;
    nodeCoor = 0.5 * sum(msh0.nodeCoor(:, [I(edgeID), J(edgeID)]), 2);
    msh.nodeCoor(:, nodeID) = nodeCoor;
    msh.nNodes = msh.nNodes + 1;
    
    edgeHasVisited(edgeID) = true;
    newEdgeNodeIDs(edgeID) = nodeID;
end
    
end

% refine an edge
function [nodeID, nodeCoor, edgeIDs] = refineEdge2(edgeID)

if edgeHasVisited(edgeID)
    nodeID   = newEdgeNodeIDs(edgeID);
    nodeCoor = msh.nodeCoor(:, nodeID);
    edgeIDs  = newEdgeIDs(:, edgeID);
else
    edgeType   = msh0.edgeType(edgeID);
    edgeNodes  = msh0.edgeNodes(:, edgeID);
    edgeNElems = msh0.edgeNElems(edgeID);
    
    nodeID = msh.nNodes + 1;
    nodeCoor = getEdgeNode(msh0.nodeCoor(:, edgeNodes(1)), msh0.nodeCoor(:, edgeNodes(2)), edgeType); 
    msh.nodeCoor(:, nodeID) = nodeCoor;    
    msh.nNodes = msh.nNodes + 1;  

    edgeIDs = [msh.nEdges + 1; msh.nEdges + 2];
    msh.edgeType(edgeIDs)     = edgeType;
    msh.edgeNodes(:, edgeIDs) = [edgeNodes(1), nodeID; nodeID, edgeNodes(2)];
    msh.edgeNElems(edgeIDs)   = edgeNElems;
    msh.nEdges = msh.nEdges + 2;
    
    edgeHasVisited(edgeID) = true;
    newEdgeNodeIDs(edgeID) = nodeID;
    newEdgeIDs(:, edgeID)  = edgeIDs;
end
   
end

% refine a face
function [nodeID, nodeCoor, faceIDs] = refineFace1(elemID, faceID, elemNodes4, elemNodeCoors4)
    
faceType  = msh0.faceType(faceID);
faceNodes = msh0.faceNodes(:, faceID);
faceElems = msh0.faceElems(:, faceID);
faceR2L   = msh0.faceR2L(faceID);
if faceHasVisited(faceID)
    nodeID   = newFaceNodeIDs(faceID);
    nodeCoor = msh.nodeCoor(:, nodeID);
    faceIDs  = newFaceIDs(:, faceID);
else
    faceEdges     = faceEdgeIDs(:, faceID);
    faceNums      = msh0.faceNums(:, faceID);
    faceNodeCoors = msh0.nodeCoor(:, faceNodes);
    % Warning : This is wrong for edges who lie on curved surface
    nodeIDs_ = zeros(1, 9);
    nodeIDs_([1, 3, 9, 7]) = faceNodes;
    nodeIDs_(4) = refineEdge1(faceEdges(1));
    nodeIDs_(6) = refineEdge1(faceEdges(2));
    nodeIDs_(2) = refineEdge1(faceEdges(3));
    nodeIDs_(8) = refineEdge1(faceEdges(4));   
    nodeIDs_(5) = msh.nNodes + 1;
    nodeID = nodeIDs_(5);
    nodeCoor = getFaceNode(faceNodeCoors(:, 1), faceNodeCoors(:, 2), faceNodeCoors(:, 3), faceNodeCoors(:, 4), faceType);
    msh.nodeCoor(:, nodeID) = nodeCoor; 
    msh.nNodes = msh.nNodes + 1;
    
    faceIDs    = msh.nFaces + 1 : msh.nFaces + 4;
    nodeCoors_ = msh.nodeCoor(:, nodeIDs_);
    msh.faceNodes(1, faceIDs) = nodeIDs_([1, 2, 4, 5]);
    msh.faceNodes(2, faceIDs) = nodeIDs_([2, 3, 5, 6]);
    msh.faceNodes(3, faceIDs) = nodeIDs_([5, 6, 8, 9]);
    msh.faceNodes(4, faceIDs) = nodeIDs_([4, 5, 7, 8]);    
    msh.faceElems(1, faceIDs) = 8 * (faceElems(1) - 1) + getLocalSubelementIDs(faceNums(1), faceR2L, true);
    msh.faceElems(2, faceIDs) = 8 * (faceElems(2) - 1) + getLocalSubelementIDs(faceNums(2), faceR2L, false);
    msh.faceNums( 1, faceIDs) = faceNums(1);
    msh.faceNums( 2, faceIDs) = faceNums(2);
    msh.faceR2L(faceIDs)      = faceR2L;
    tempNormals = 0.5 * mycross3(nodeCoors_(:, [1, 2, 4, 5]) - nodeCoors_(:, [5, 6, 8, 9]), nodeCoors_(:,[2, 3, 5, 6]) - nodeCoors_(:, [4, 5, 7, 8]));
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    if mod(faceNums(1), 2) == 1
        tempNormals = -tempNormals;
    end       
    msh.faceNormalx(faceIDs) = tempNormals(1, :);
    msh.faceNormaly(faceIDs) = tempNormals(2, :);
    msh.faceNormalz(faceIDs) = tempNormals(3, :);
    msh.faceSize(faceIDs)    = tempSizes;
    msh.faceType(faceIDs)    = faceType;
    msh.faceJac(faceIDs)     = tempSizes / 4;    
    if msh0.faceIsCurved(faceID)
        msh.faceIsCurved(faceIDs) = true;
        msh.elemIsCurved(msh.faceElems(1, faceIDs)) = true;
    end
    msh.nFaces = msh.nFaces + 4;      
    
    faceHasVisited(faceID) = true;
    newFaceNodeIDs(faceID) = nodeID;
    newFaceIDs(:, faceID)  = faceIDs;
end
if (faceElems(2) == elemID)
    faceIDs = reorderFaceIDs_l2r(faceIDs, faceR2L);
end
if (faceType == 1)
    if ~(all(faceNodes == elemNodes4) || all(faceNodes == reorderNodeIDs_r2l(elemNodes4, faceR2L)))
        nodeID = msh.nNodes + 1;
        nodeCoor = 0.25 * sum(elemNodeCoors4, 2);  
        msh.nodeCoor(:, nodeID) = nodeCoor;        
        msh.nNodes = msh.nNodes + 1;         
    end  
end

end

% refine a face
function [nodeID, nodeCoor, edgeIDs, faceIDs] = refineFace2(elemID, faceID, elemNodes9, elemNodeCoors4)
    
faceType  = msh0.faceType(faceID);
faceNodes = msh0.faceNodes(:, faceID);
faceElems = msh0.faceElems(:, faceID);
faceR2L   = msh0.faceR2L(faceID);    
if faceHasVisited(faceID)
    nodeID   = newFaceNodeIDs(faceID);
    nodeCoor = msh.nodeCoor(:, nodeID);
    edgeIDs  = newFaceEdgeIDs(:, faceID);
    faceIDs  = newFaceIDs(:, faceID);
else
    faceEdges     = msh0.faceEdges(:, faceID);
    faceNums      = msh0.faceNums(:, faceID);
    faceNodeCoors = msh0.nodeCoor(:, faceNodes);
    faceEdgeTypes = msh0.edgeType(faceEdges);
    faceEdgeNodes = msh0.edgeNodes(:, faceEdges);
    
    nodeIDs_  = zeros(1, 9);
    edgeIDs1_ = zeros(3, 2);
    edgeIDs2_ = zeros(2, 3);
   
    nodeIDs_([1, 3, 9, 7]) = faceNodes;
    [nodeIDs_(4), ~, edgeIDs1_(1, :)] = refineEdge2(faceEdges(1));
    [nodeIDs_(6), ~, edgeIDs1_(3, :)] = refineEdge2(faceEdges(2));
    [nodeIDs_(2), ~, edgeIDs2_(:, 1)] = refineEdge2(faceEdges(3));
    [nodeIDs_(8), ~, edgeIDs2_(:, 3)] = refineEdge2(faceEdges(4));
    if (faceEdgeNodes(1, 1) == nodeIDs_(7))
        edgeIDs1_(1, :) = flip(edgeIDs1_(1, :));
    end
    if (faceEdgeNodes(1, 2) == nodeIDs_(9))
        edgeIDs1_(3, :) = flip(edgeIDs1_(3, :));
    end   
    if (faceEdgeNodes(1, 3) == nodeIDs_(3))
        edgeIDs2_(:, 1) = flip(edgeIDs2_(:, 1));
    end  
    if (faceEdgeNodes(1, 4) == nodeIDs_(9))
        edgeIDs2_(:, 3) = flip(edgeIDs2_(:, 3));
    end    
    
    nodeIDs_(5) = msh.nNodes + 1;
    nodeID = nodeIDs_(5);
    nodeCoor = getFaceNode(faceNodeCoors(:, 1), faceNodeCoors(:, 2), faceNodeCoors(:, 3), faceNodeCoors(:, 4), faceType);
    msh.nodeCoor(:, nodeIDs_(5)) = nodeCoor;
    msh.nNodes = msh.nNodes + 1;

    edgeIDs = msh.nEdges + 1 : msh.nEdges + 4;
    edgeIDs1_(2, :) = edgeIDs(1 : 2);
    edgeIDs2_(:, 2) = edgeIDs(3 : 4)';
    if all(faceEdgeTypes(1) == faceEdgeTypes(2 : 4))
        msh.edgeType(edgeIDs) = faceEdgeTypes(1);
    end    
    msh.edgeNodes(:, edgeIDs) = [nodeIDs_(2), nodeIDs_(5), nodeIDs_(4), nodeIDs_(5); ...
                                 nodeIDs_(5), nodeIDs_(8), nodeIDs_(5), nodeIDs_(6)];                
    if (faceType == 0)
        msh.edgeNElems(edgeIDs) = 4;
    else
        msh.edgeNElems(edgeIDs) = 2;
    end
    msh.nEdges = msh.nEdges + 4;    
    
    faceIDs    = msh.nFaces + 1 : msh.nFaces + 4;
    nodeCoors_ = msh.nodeCoor(:, nodeIDs_);
    msh.faceNodes(1, faceIDs) = nodeIDs_([1, 2, 4, 5]);
    msh.faceNodes(2, faceIDs) = nodeIDs_([2, 3, 5, 6]);
    msh.faceNodes(3, faceIDs) = nodeIDs_([5, 6, 8, 9]);
    msh.faceNodes(4, faceIDs) = nodeIDs_([4, 5, 7, 8]);  
    msh.faceEdges(1, faceIDs) = reshape(edgeIDs1_(1 : 2, :), [1, 4]);
    msh.faceEdges(2, faceIDs) = reshape(edgeIDs1_(2 : 3, :), [1, 4]);
    msh.faceEdges(3, faceIDs) = reshape(edgeIDs2_(:, 1 : 2), [1, 4]);
    msh.faceEdges(4, faceIDs) = reshape(edgeIDs2_(:, 2 : 3), [1, 4]);    
    msh.faceElems(1, faceIDs) = 8 * (faceElems(1) - 1) + getLocalSubelementIDs(faceNums(1), faceR2L, true);
    msh.faceElems(2, faceIDs) = 8 * (faceElems(2) - 1) + getLocalSubelementIDs(faceNums(2), faceR2L, false);
    msh.faceNums( 1, faceIDs) = faceNums(1);
    msh.faceNums( 2, faceIDs) = faceNums(2);
    msh.faceR2L(faceIDs)      = faceR2L;
    tempNormals = 0.5 * mycross3(nodeCoors_(:, [1, 2, 4, 5]) - nodeCoors_(:, [5, 6, 8, 9]), nodeCoors_(:,[2, 3, 5, 6]) - nodeCoors_(:, [4, 5, 7, 8]));
    tempSizes   = vecnorm(tempNormals);
    tempNormals = tempNormals ./ tempSizes;
    if mod(faceNums(1), 2) == 1
        tempNormals = -tempNormals;
    end       
    msh.faceNormalx(faceIDs) = tempNormals(1, :);
    msh.faceNormaly(faceIDs) = tempNormals(2, :);
    msh.faceNormalz(faceIDs) = tempNormals(3, :);
    msh.faceSize(faceIDs)    = tempSizes;
    msh.faceType(faceIDs)    = faceType;
    msh.faceJac(faceIDs)     = tempSizes / 4;    
    if msh0.faceIsCurved(faceID)
        msh.faceIsCurved(faceIDs) = true;
        msh.elemIsCurved(msh.faceElems(1, faceIDs)) = true;
    end
    msh.nFaces = msh.nFaces + 4;  
    
    faceHasVisited(faceID) = true;
    newFaceNodeIDs(faceID) = nodeID;
    newFaceEdgeIDs(:, faceID) = edgeIDs;
    newFaceIDs(:, faceID) = faceIDs;    
end
if (faceElems(2) == elemID)
    edgeIDs = reorderEdgeIDs_l2r(edgeIDs, faceR2L);
    faceIDs = reorderFaceIDs_l2r(faceIDs, faceR2L); 
end
if (faceType == 1)
    if ~(all(faceNodes == elemNodes9([1, 3, 9, 7])) || all(faceNodes == reorderNodeIDs_r2l(elemNodes9([1, 3, 9, 7]), faceR2L)))
        nodeID = msh.nNodes + 1;
        nodeCoor = 0.25 * sum(elemNodeCoors4, 2);  
        msh.nodeCoor(:, nodeID) = nodeCoor;        
        msh.nNodes = msh.nNodes + 1;   
        
        edgeIDs = msh.nEdges + 1 : msh.nEdges + 4;
        msh.edgeType(edgeIDs)     = 0;
        msh.edgeNodes(:, edgeIDs) = [elemNodes9(2),        nodeID, elemNodes9(4),         nodeID; ...
                                            nodeID, elemNodes9(8),        nodeID, elemNodes9(6)];  
        msh.edgeNElems(edgeIDs)   = 2; 
        msh.nEdges = msh.nEdges + 4;
    end  
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

% Given the four face node IDs ordered w.r.t. the rifght element, and the
% face-to-face index, reorder them according to the left element
function IDs = reorderNodeIDs_r2l(IDs, f2f)
    
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

% Given the internal face edge IDs ordered w.r.t. the left element, and the
% face-to-face index, reorder them according to the right element
function IDs = reorderEdgeIDs_l2r(IDs, f2f)
    
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

function subeIDs = getLocalSubelementIDs(fn, f2f, isLeft)
    
switch fn
    case 1
        subeIDs = [1, 3, 5, 7];
    case 2
        subeIDs = [2, 4, 6, 8];
    case 3
        subeIDs = [1, 5, 2, 6];
    case 4
        subeIDs = [3, 7, 4, 8];
    case 5
        subeIDs = 1 : 4;
    case 6
        subeIDs = 5 : 8;
    otherwise
        error('Wrong face ID')
end
if ~isLeft
    subeIDs = reorderFaceIDs_l2r(subeIDs, f2f);
end

end





