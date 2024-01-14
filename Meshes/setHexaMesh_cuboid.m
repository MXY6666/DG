% The ordering of the vertices and faces in a cube:
%
%     5o---------o8
%     /|        /|
%    / |       / |
%  6o---------o7 |
%   | 1o------|--o4
%   | /       | /
%   |/        |/
%  2o---------o3
%
% Face1 : 1485; Face2 : 2376; 
% Face3 : 1562; Face4 : 4873;
% Face5 : 1234; Face6 : 5678;
% Edge1 : 12; Edge2  : 23; Edge3  : 34; Edge4  : 41;
% Edge5 : 15; Edge6  : 26; Edge7  : 37; Edge8  : 48;
% Edge9 : 56; Edge10 : 67; Edge11 : 78; Edge12 : 85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain     : a 6D vector describing the computational domain 
% N          : 3D vector, number of elements in the x, y and z direction
% bcs        : 6D vector containing the sixth boundary types
% type = 311 : ignore the adaptive information      
% type = 312 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% quad       : Gauss quadrature rule for a hexahedron
% bs         : basis function set
% Copyright (C) Xiangyi Meng.

function msh = setHexaMesh_cuboid(domain, N, bcs, type, maxLevel, quad, bs)

if (nargin < 2)
    error('Not enough arguments')
end

if (length(domain) ~= 6)
    error('Wrong size of argument domain')
end

if (length(N) ~= 3)
    error('Wrong size of argument N')
end

if (domain(2) <= domain(1)) || (domain(4) <= domain(3)) || (domain(6) <= domain(5))
    error('Wrong argument domain')
end

if (nargin < 3)  || isempty(bcs)
    bcs = ones(1, 6);
end
if (length(bcs) ~= 6)
    error('Wrong boundary conditions')
end
if (bcs(1) == 1)
    if (bcs(2) == 1)
        isPeriodicInX = true;
    else
        error('Wrong boundary condition in x direction')
    end
else
     if (bcs(2) == 1)
         error('Wrong boundary condition in x direction')
     else
         isPeriodicInX = false;
     end     
end
if (bcs(3) == 1)
    if (bcs(4) == 1)
        isPeriodicInY = true;
    else
        error('Wrong boundary condition in y direction')
    end
else
     if (bcs(4) == 1)
         error('Wrong boundary condition in y direction')
     else
         isPeriodicInY = false;
     end   
end
if (bcs(5) == 1)
    if (bcs(6) == 1)
        isPeriodicInZ = true;
    else
        error('Wrong boundary condition in z direction')
    end
else
     if (bcs(6) == 1)
         error('Wrong boundary condition in z direction')
     else
         isPeriodicInZ = false;
     end   
end

if (nargin < 4)  || isempty(type)
    type = 311;
end
if (type ~= 311) && (type ~= 312) 
    error('wrong mesh type')
end

if (nargin < 5)  || isempty(maxLevel) || (type == 311) 
    maxLevel = 0;
end
if (type == 312) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 6) || isempty(quad)
    quad = GaussQuadratureRule_hexa([2, 2, 2], 301);
end

if (nargin < 7) || isempty(bs)
    quad2 = GaussQuadratureRule_square([2, 2], 201);
    bs = setBasisFunctionSet_cube(quad2, quad, 301);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = prod(N);
nn = prod(N + 1);
n12 = N(1) * N(2); n23 = N(2) * N(3); n31 = N(3) * N(1); 
ind_e = reshape(1 : ne, N);
ind_f1 = zeros(N(1) + 1, N(2), N(3));
ind_f2 = zeros(N(1), N(2) + 1, N(3));
ind_f3 = zeros(N(1), N(2), N(3) + 1);
if isPeriodicInX
    ind_f1(1 : end - 1, :, :) = reshape(1 : ne, N);
    ind_f1(end, :, :) = ind_f1(1, :, :);
    nf = ne;
else
    ind_f1 = reshape(1 : ne + n23, size(ind_f1));
    nf = ne + n23;
end
if isPeriodicInY
    ind_f2(:, 1 : end - 1, :) = reshape(nf + 1 : nf + ne, N);
    ind_f2(:, end, :) = ind_f2(:, 1, :);
    nf = nf + ne;
else
    ind_f2 = reshape(nf + 1 : nf + ne + n31, size(ind_f2));
    nf = nf + ne + n31;
end
if isPeriodicInZ
    ind_f3(:, :, 1 : end - 1) = reshape(nf + 1 : nf + ne, N);
    ind_f3(:, :, end) = ind_f3(:, :, 1);
else
    ind_f3 = reshape(nf + 1 : nf + ne + n12, size(ind_f3));
end
if (type == 312)    
    ind_ed1 = reshape(1 : ne + n12 + n31 + N(1), [N(1), N(2) + 1, N(3) + 1]);
    
    ned = ne + n12 + n31 + N(1);
    ind_ed2 = reshape(ned + 1 : ned + ne + n12 + n23 + N(2), [N(1) + 1, N(2), N(3) + 1]);
    
    ned = ned + ne + n12 + n23 + N(2);
    ind_ed3 = reshape(ned + 1 : ned + ne + n31 + n23 + N(3), [N(1) + 1, N(2) + 1, N(3)]);
end
ind_n = reshape(1 : nn, N + 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msh          = getEmptyMesh;
msh.dm       = domain;
msh.N        = N;
msh.type     = type;
msh.maxLevel = maxLevel;
msh.bndTypes = unique(bcs, 'stable');
msh.nElems   = ne;
msh.nFaces   = 3 * ne + n23 + n31 + n12; 
if isPeriodicInX
    msh.nFaces = msh.nFaces - n23;
end
if isPeriodicInY
    msh.nFaces = msh.nFaces - n31;
end
if isPeriodicInZ
    msh.nFaces = msh.nFaces - n12;
end
if (type == 312)
   msh.nEdges = 3 * ne + 2 * (n12 + n23 + n31) + sum(N); 
end
msh.nNodes  = nn;
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems  = 1 : msh.nElems;
msh.LFaces  = 1 : msh.nFaces;
msh.nGPs    = quad.nps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of nodes
msh.nodeCoor = zeros(3, msh.nNodes);
msh.nodeCoor(1, :) = repmat(linspace(domain(1), domain(2), N(1) + 1), [1, (N(2) + 1) * (N(3) + 1)]);
msh.nodeCoor(2, :) = repmat(repelem(linspace(domain(3), domain(4), N(2) + 1), N(1) + 1), [1, N(3) + 1]);
msh.nodeCoor(3, :) = repelem(linspace(domain(5), domain(6), N(3) + 1), (N(1) + 1) * (N(2) + 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of elements
% (some) topology information of elements
msh.elemNodes = zeros(8, msh.nElems);
msh.elemFaces = zeros(6, msh.nElems);

msh.elemNodes(1, :) = reshape(ind_n(1 : end - 1, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(2, :) = reshape(ind_n(2 : end, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(3, :) = reshape(ind_n(2 : end, 2 : end, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(4, :) = reshape(ind_n(1 : end - 1, 2 : end, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(5, :) = reshape(ind_n(1 : end - 1, 1 : end - 1, 2 : end), [1, msh.nElems]);
msh.elemNodes(6, :) = reshape(ind_n(2 : end, 1 : end - 1, 2 : end), [1, msh.nElems]);
msh.elemNodes(7, :) = reshape(ind_n(2 : end, 2 : end, 2 : end), [1, msh.nElems]);
msh.elemNodes(8, :) = reshape(ind_n(1 : end - 1, 2 : end, 2 : end), [1, msh.nElems]);
msh.elemFaces(1, :) = reshape(ind_f1(1 : end - 1, :, :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(2 : end, :, :), [1, msh.nElems]);
msh.elemFaces(3, :) = reshape(ind_f2(:, 1 : end - 1, :), [1, msh.nElems]);
msh.elemFaces(4, :) = reshape(ind_f2(:, 2 : end, :), [1, msh.nElems]);
msh.elemFaces(5, :) = reshape(ind_f3(:, :, 1 : end - 1), [1, msh.nElems]);
msh.elemFaces(6, :) = reshape(ind_f3(:, :, 2 : end), [1, msh.nElems]);

% (some) metric information of elements
x = reshape(msh.nodeCoor(1, msh.elemNodes), [8, msh.nElems]);
y = reshape(msh.nodeCoor(2, msh.elemNodes), [8, msh.nElems]);
z = reshape(msh.nodeCoor(3, msh.elemNodes), [8, msh.nElems]);
msh.elemSize = 1 / 6 * dot(msh.nodeCoor(:, msh.elemNodes(7, :)) - msh.nodeCoor(:, msh.elemNodes(1, :)), mycross3(msh.nodeCoor(:, msh.elemNodes(2, :)) - msh.nodeCoor(:, msh.elemNodes(1, :)), msh.nodeCoor(:, msh.elemNodes(3, :)) - msh.nodeCoor(:, msh.elemNodes(6, :))) + mycross3(msh.nodeCoor(:, msh.elemNodes(4, :)) - msh.nodeCoor(:, msh.elemNodes(1, :)), msh.nodeCoor(:, msh.elemNodes(8, :)) - msh.nodeCoor(:, msh.elemNodes(3, :))) + mycross3(msh.nodeCoor(:, msh.elemNodes(5, :)) - msh.nodeCoor(:, msh.elemNodes(1, :)), msh.nodeCoor(:, msh.elemNodes(6, :)) - msh.nodeCoor(:, msh.elemNodes(8, :))));
msh.elemDiam = msh.elemSize.^(1/3);

% other useful information for computation
xxi   = quad.N_xi * x;
yxi   = quad.N_xi * y;
zxi   = quad.N_xi * z;
xeta  = quad.N_eta * x;
yeta  = quad.N_eta * y;
zeta  = quad.N_eta * z;
xzeta = quad.N_zeta * x;
yzeta = quad.N_zeta * y;
zzeta = quad.N_zeta * z;

msh.elemGPx = quad.N * x;
msh.elemGPy = quad.N * y;
msh.elemGPz = quad.N * z;
msh.elemJac = xxi   .* (yeta  .* zzeta - yzeta .* zeta) ...
            + xeta  .* (yzeta .* zxi   - yxi   .* zzeta) ...
            + xzeta .* (yxi   .* zeta  - yeta  .* zxi);
if any(msh.elemJac < 0)
    error('Negative Jacobian')
end
msh.elemJxix   = yeta  .* zzeta - yzeta .* zeta;
msh.elemJxiy   = zeta  .* xzeta - zzeta .* xeta;
msh.elemJxiz   = xeta  .* yzeta - xzeta .* yeta;
msh.elemJetax  = yzeta .* zxi   - yxi   .* zzeta;
msh.elemJetay  = zzeta .* xxi   - zxi   .* xzeta;
msh.elemJetaz  = xzeta .* yxi   - xxi   .* yzeta;
msh.elemJzetax = yxi   .* zeta  - yeta  .* zxi;
msh.elemJzetay = zxi   .* xeta  - zeta  .* xxi;
msh.elemJzetaz = xxi   .* yeta  - xeta  .* yxi;
msh.elemIsCurved = false(1, msh.nElems);

% Store inverse of mass matrix in each element
massMat        = bs.phitw{1} * (repmat(bs.phi{1}, [1, msh.nElems]) .* repelem(msh.elemJac, 1, bs.nb));
massMatInv     = cellfun(@inv, mat2cell(massMat, bs.nb, bs.nb * ones(1, msh.nElems)), 'UniformOutput', false);
msh.massMatInv = cellfun(@mychop, massMatInv, num2cell(1.0e-12 * ones(1, msh.nElems)), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of faces
% (some) topology information of faces
msh.faceNodes = zeros(4, msh.nFaces);
msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums  = zeros(2, msh.nFaces); 
msh.faceR2L   = ones(1, msh.nFaces);

msh.faceNodes(1,  ind_f1(1 : end - 1, :, :)) = reshape(ind_n(1 : end - 1, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(2,  ind_f1(1 : end - 1, :, :)) = reshape(ind_n(1 : end - 1, 2 : end, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(3,  ind_f1(1 : end - 1, :, :)) = reshape(ind_n(1 : end - 1, 2 : end, 2 : end), [1, msh.nElems]);
msh.faceNodes(4,  ind_f1(1 : end - 1, :, :)) = reshape(ind_n(1 : end - 1, 1 : end - 1, 2 : end), [1, msh.nElems]);
msh.faceElems(1,  ind_f1(2 : end - 1, :, :)) = reshape(ind_e(1 : end - 1, :, :), [1, msh.nElems - n23]);
msh.faceElems(2,  ind_f1(2 : end - 1, :, :)) = reshape(ind_e(2 : end, :, :), [1, msh.nElems - n23]);
msh.faceNums( 1,  ind_f1(2 : end - 1, :, :)) = 2;
msh.faceNums( 2,  ind_f1(2 : end - 1, :, :)) = 1;

msh.faceNodes(1,  ind_f2(:, 1 : end - 1, :)) = reshape(ind_n(1 : end - 1, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(2,  ind_f2(:, 1 : end - 1, :)) = reshape(ind_n(1 : end - 1, 1 : end - 1, 2 : end), [1, msh.nElems]);
msh.faceNodes(3,  ind_f2(:, 1 : end - 1, :)) = reshape(ind_n(2 : end, 1 : end - 1, 2 : end), [1, msh.nElems]);
msh.faceNodes(4,  ind_f2(:, 1 : end - 1, :)) = reshape(ind_n(2 : end, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceElems(1,  ind_f2(:, 2 : end - 1, :)) = reshape(ind_e(:, 1 : end - 1, :), [1, msh.nElems - n31]);
msh.faceElems(2,  ind_f2(:, 2 : end - 1, :)) = reshape(ind_e(:, 2 : end, :), [1, msh.nElems - n31]);
msh.faceNums( 1,  ind_f2(:, 2 : end - 1, :)) = 4;
msh.faceNums( 2,  ind_f2(:, 2 : end - 1, :)) = 3;

msh.faceNodes(1,  ind_f3(:, :, 1 : end - 1)) = reshape(ind_n(1 : end - 1, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(2,  ind_f3(:, :, 1 : end - 1)) = reshape(ind_n(2 : end, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(3,  ind_f3(:, :, 1 : end - 1)) = reshape(ind_n(2 : end, 2 : end, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(4,  ind_f3(:, :, 1 : end - 1)) = reshape(ind_n(1 : end - 1, 2 : end, 1 : end - 1), [1, msh.nElems]);
msh.faceElems(1,  ind_f3(:, :, 2 : end - 1)) = reshape(ind_e(:, :, 1 : end - 1), [1, msh.nElems - n12]);
msh.faceElems(2,  ind_f3(:, :, 2 : end - 1)) = reshape(ind_e(:, :, 2 : end), [1, msh.nElems - n12]);
msh.faceNums( 1,  ind_f3(:, :, 2 : end - 1)) = 6;
msh.faceNums( 2,  ind_f3(:, :, 2 : end - 1)) = 5;

if isPeriodicInX 
    msh.faceElems(1, ind_f1(1, :, :)) = reshape(ind_e(end, :, :), [1, n23]);
    msh.faceElems(2, ind_f1(1, :, :)) = reshape(ind_e(1, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(1, :, :)) = 2;
    msh.faceNums( 2, ind_f1(1, :, :)) = 1;    
else
    msh.faceElems(1, ind_f1(1, :, :)) = reshape(ind_e(1, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(1, :, :)) = 1;
    msh.faceNodes(1, ind_f1(end, :, :)) = reshape(ind_n(end, 1 : end - 1, 1 : end - 1), [1, n23]);
    msh.faceNodes(2, ind_f1(end, :, :)) = reshape(ind_n(end, 2 : end, 1 : end - 1), [1, n23]);
    msh.faceNodes(3, ind_f1(end, :, :)) = reshape(ind_n(end, 2 : end, 2 : end), [1, n23]);
    msh.faceNodes(4, ind_f1(end, :, :)) = reshape(ind_n(end, 1 : end - 1, 2 : end), [1, n23]);     
    msh.faceElems(1, ind_f1(end, :, :)) = reshape(ind_e(end, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(end, :, :)) = 2;
end
if isPeriodicInY
    msh.faceElems(1, ind_f2(:, 1, :)) = reshape(ind_e(:, end, :), [1, n31]);
    msh.faceElems(2, ind_f2(:, 1, :)) = reshape(ind_e(:, 1, :), [1, n31]);
    msh.faceNums( 1, ind_f2(:, 1, :)) = 4;
    msh.faceNums( 2, ind_f2(:, 1, :)) = 3;
else
    msh.faceElems(1, ind_f2(:, 1, :)) = reshape(ind_e(:, 1, :), [1, n31]);
    msh.faceNums( 1, ind_f2(:, 1, :)) = 3;
    msh.faceNodes(1, ind_f2(:, end, :)) = reshape(ind_n(1 : end - 1, end, 1 : end - 1), [1, n31]);
    msh.faceNodes(2, ind_f2(:, end, :)) = reshape(ind_n(1 : end - 1, end, 2 : end), [1, n31]);
    msh.faceNodes(3, ind_f2(:, end, :)) = reshape(ind_n(2 : end, end, 2 : end), [1, n31]);
    msh.faceNodes(4, ind_f2(:, end, :)) = reshape(ind_n(2 : end, end, 1 : end - 1), [1, n31]);    
    msh.faceElems(1, ind_f2(:, end, :)) = reshape(ind_e(:, end, :), [1, n31]);
    msh.faceNums( 1, ind_f2(:, end, :)) = 4;    
end
if isPeriodicInZ
    msh.faceElems(1, ind_f3(:, :, 1)) = reshape(ind_e(:, :, end), [1, n12]);
    msh.faceElems(2, ind_f3(:, :, 1)) = reshape(ind_e(:, :, 1), [1, n12]);
    msh.faceNums( 1, ind_f3(:, :, 1)) = 6;
    msh.faceNums( 2, ind_f3(:, :, 1)) = 5;
else
    msh.faceElems(1, ind_f3(:, :, 1)) = reshape(ind_e(:, :, 1), [1, n12]);
    msh.faceNums( 1, ind_f3(:, :, 1)) = 5;
    msh.faceNodes(1, ind_f3(:, :, end)) = reshape(ind_n(1 : end - 1, 1 : end - 1, end), [1, n12]);
    msh.faceNodes(2, ind_f3(:, :, end)) = reshape(ind_n(2 : end, 1 : end - 1, end), [1, n12]);
    msh.faceNodes(3, ind_f3(:, :, end)) = reshape(ind_n(2 : end, 2 : end, end), [1, n12]);
    msh.faceNodes(4, ind_f3(:, :, end)) = reshape(ind_n(1 : end - 1, 2 : end, end), [1, n12]);      
    msh.faceElems(1, ind_f3(:, :, end)) = reshape(ind_e(:, :, end), [1, n12]);
    msh.faceNums( 1, ind_f3(:, :, end)) = 6;  
end

% (some) metric information of faces
ds = 0.5 * mycross3(msh.nodeCoor(:, msh.faceNodes(1, :)) - msh.nodeCoor(:, msh.faceNodes(3, :)), msh.nodeCoor(:, msh.faceNodes(2, :)) - msh.nodeCoor(:, msh.faceNodes(4, :)));
msh.faceSize = vecnorm(ds);
msh.faceNormalx = ds(1, :) ./ msh.faceSize;
msh.faceNormaly = ds(2, :) ./ msh.faceSize;
msh.faceNormalz = ds(3, :) ./ msh.faceSize;
msh.faceType = zeros(1, msh.nFaces);
msh.faceType(ind_f1(1, :, :)) = bcs(1);
msh.faceType(ind_f2(:, 1, :)) = bcs(3);
msh.faceType(ind_f3(:, :, 1)) = bcs(5);
if ~isPeriodicInX 
    msh.faceNormalx(ind_f1(1, :, :)) = -msh.faceNormalx(ind_f1(1, :, :));
    msh.faceType(ind_f1(end, :, :)) = bcs(2);
end
if ~isPeriodicInY
    msh.faceNormaly(ind_f2(:, 1, :)) = -msh.faceNormaly(ind_f2(:, 1, :));
    msh.faceType(ind_f2(:, end, :)) = bcs(4); 
end
if ~isPeriodicInZ
    msh.faceNormalz(ind_f3(:, :, 1)) = -msh.faceNormalz(ind_f3(:, :, 1));
    msh.faceType(ind_f3(:, :, end)) = bcs(6); 
end

% other useful information for computation
msh.faceJac      = 0.25 * msh.faceSize;
msh.faceIsCurved = false(1, msh.nFaces);

% Classify internal faces and boundary faces
msh.intLFaces  = find(msh.faceType == 0);
msh.nIntLFaces = length(msh.intLFaces);

msh.bndLFaces  = cell(1, length(msh.bndTypes));
msh.nBndLFaces = zeros(1, length(msh.bndTypes));
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = find(msh.faceType == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (type == 312)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the fields of edges
    % metric information of edges
    msh.edgeType = zeros(1, msh.nEdges);
    
    % topology information of edges
    msh.edgeNodes(1, ind_ed1) = reshape(ind_n(1 : end - 1, :, :), [1, msh.nElems + n12 + n31 + N(1)]);
    msh.edgeNodes(2, ind_ed1) = reshape(ind_n(2 : end, :, :), [1, msh.nElems + n12 + n31 + N(1)]);
    msh.edgeNodes(1, ind_ed2) = reshape(ind_n(:, 1 : end - 1, :), [1, msh.nElems + n12 + n23 + N(2)]);
    msh.edgeNodes(2, ind_ed2) = reshape(ind_n(:, 2 : end, :), [1, msh.nElems + n12 + n23 + N(2)]);
    msh.edgeNodes(1, ind_ed3) = reshape(ind_n(:, :, 1 : end - 1), [1, msh.nElems + n31 + n23 + N(3)]);
    msh.edgeNodes(2, ind_ed3) = reshape(ind_n(:, :, 2 : end), [1, msh.nElems + n31 + n23 + N(3)]);
    
    msh.edgeNElems(ind_ed1(:, 1, 1))                     = 1;
    msh.edgeNElems(ind_ed1(:, 1, end))                   = 1;
    msh.edgeNElems(ind_ed1(:, end, 1))                   = 1;  
    msh.edgeNElems(ind_ed1(:, end, end))                 = 1;
    msh.edgeNElems(ind_ed1(:, 1, 2 : end - 1))           = 2;
    msh.edgeNElems(ind_ed1(:, end, 2 : end - 1))         = 2;
    msh.edgeNElems(ind_ed1(:, 2 : end - 1, 1))           = 2;
    msh.edgeNElems(ind_ed1(:, 2 : end - 1, end))         = 2;   
    msh.edgeNElems(ind_ed1(:, 2 : end - 1, 2 : end - 1)) = 4; 
    msh.edgeNElems(ind_ed2(1, :, 1))                     = 1;
    msh.edgeNElems(ind_ed2(1, :, end))                   = 1;
    msh.edgeNElems(ind_ed2(end, :, 1))                   = 1;
    msh.edgeNElems(ind_ed2(end, :, end))                 = 1;
    msh.edgeNElems(ind_ed2(1, :, 2 : end - 1))           = 2;
    msh.edgeNElems(ind_ed2(end, :, 2 : end - 1))         = 2;
    msh.edgeNElems(ind_ed2(2 : end - 1, :, 1))           = 2;
    msh.edgeNElems(ind_ed2(2 : end - 1, :, end))         = 2;
    msh.edgeNElems(ind_ed2(2 : end - 1, :, 2 : end - 1)) = 4;
    msh.edgeNElems(ind_ed3(1, 1, :))                     = 1;
    msh.edgeNElems(ind_ed3(1, end, :))                   = 1;
    msh.edgeNElems(ind_ed3(end, 1, :))                   = 1;
    msh.edgeNElems(ind_ed3(end, end, :))                 = 1;
    msh.edgeNElems(ind_ed3(1, 2 : end - 1, :))           = 2;
    msh.edgeNElems(ind_ed3(end, 2 : end - 1, :))         = 2;
    msh.edgeNElems(ind_ed3(2 : end - 1, 1, :))           = 2;
    msh.edgeNElems(ind_ed3(2 : end - 1, end, :))         = 2;
    msh.edgeNElems(ind_ed3(2 : end - 1, 2 : end - 1, :)) = 4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Augment the fields of elements with egde information
    msh.elemEdges(1, :)  = reshape(ind_ed1(:, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
    msh.elemEdges(2, :)  = reshape(ind_ed1(:, 2 : end, 1 : end - 1), [1, msh.nElems]);
    msh.elemEdges(3, :)  = reshape(ind_ed1(:, 1 : end - 1, 2 : end), [1, msh.nElems]);
    msh.elemEdges(4, :)  = reshape(ind_ed1(:, 2 : end, 2 : end), [1, msh.nElems]);
    msh.elemEdges(5, :)  = reshape(ind_ed2(1 : end - 1, :, 1 : end - 1), [1, msh.nElems]);
    msh.elemEdges(6, :)  = reshape(ind_ed2(2 : end, :, 1 : end - 1), [1, msh.nElems]);
    msh.elemEdges(7, :)  = reshape(ind_ed2(1 : end - 1, :, 2 : end), [1, msh.nElems]);
    msh.elemEdges(8, :)  = reshape(ind_ed2(2 : end, :, 2 : end), [1, msh.nElems]);
    msh.elemEdges(9, :)  = reshape(ind_ed3(1 : end - 1, 1 : end - 1, :), [1, msh.nElems]);
    msh.elemEdges(10, :) = reshape(ind_ed3(2 : end, 1 : end - 1, :), [1, msh.nElems]);
    msh.elemEdges(11, :) = reshape(ind_ed3(1 : end - 1, 2 : end, :), [1, msh.nElems]);
    msh.elemEdges(12, :) = reshape(ind_ed3(2 : end, 2 : end, :), [1, msh.nElems]);    
    
    % Augment the fields of faces with egde information
    msh.faceEdges(1, ind_f1(1 : end - 1, :, :)) = reshape(ind_ed3(1 : end - 1, 1 : end - 1, :), [1, msh.nElems]);
    msh.faceEdges(2, ind_f1(1 : end - 1, :, :)) = reshape(ind_ed3(1 : end - 1, 2 : end, :), [1, msh.nElems]);
    msh.faceEdges(3, ind_f1(1 : end - 1, :, :)) = reshape(ind_ed2(1 : end - 1, :, 1 : end - 1), [1, msh.nElems]);
    msh.faceEdges(4, ind_f1(1 : end - 1, :, :)) = reshape(ind_ed2(1 : end - 1, :, 2 : end), [1, msh.nElems]);
          
    msh.faceEdges(1, ind_f2(:, 1 : end - 1, :)) = reshape(ind_ed1(:, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
    msh.faceEdges(2, ind_f2(:, 1 : end - 1, :)) = reshape(ind_ed1(:, 1 : end - 1, 2 : end), [1, msh.nElems]);
    msh.faceEdges(3, ind_f2(:, 1 : end - 1, :)) = reshape(ind_ed3(1 : end - 1, 1 : end - 1, :), [1, msh.nElems]);
    msh.faceEdges(4, ind_f2(:, 1 : end - 1, :)) = reshape(ind_ed3(2 : end, 1 : end - 1, :), [1, msh.nElems]);
           
    msh.faceEdges(1, ind_f3(:, :, 1 : end - 1)) = reshape(ind_ed2(1 : end - 1, :, 1 : end - 1), [1, msh.nElems]);
    msh.faceEdges(2, ind_f3(:, :, 1 : end - 1)) = reshape(ind_ed2(2 : end, :, 1 : end - 1), [1, msh.nElems]);
    msh.faceEdges(3, ind_f3(:, :, 1 : end - 1)) = reshape(ind_ed1(:, 1 : end - 1, 1 : end - 1), [1, msh.nElems]);
    msh.faceEdges(4, ind_f3(:, :, 1 : end - 1)) = reshape(ind_ed1(:, 2 : end, 1 : end - 1), [1, msh.nElems]);  
    
    if ~isPeriodicInX        
        msh.faceEdges(1, ind_f1(end, :, :)) = reshape(ind_ed3(end, 1 : end - 1, :), [1, n23]); 
        msh.faceEdges(2, ind_f1(end, :, :)) = reshape(ind_ed3(end, 2 : end, :), [1, n23]);
        msh.faceEdges(3, ind_f1(end, :, :)) = reshape(ind_ed2(end, :, 1 : end - 1), [1, n23]);
        msh.faceEdges(4, ind_f1(end, :, :)) = reshape(ind_ed2(end, :, 2 : end), [1, n23]);
    end
    if ~isPeriodicInY          
        msh.faceEdges(1, ind_f2(:, end, :)) = reshape(ind_ed1(:, end, 1 : end - 1), [1, n31]); 
        msh.faceEdges(2, ind_f2(:, end, :)) = reshape(ind_ed1(:, end, 2 : end), [1, n31]);
        msh.faceEdges(3, ind_f2(:, end, :)) = reshape(ind_ed3(1 : end - 1, end, :), [1, n31]);
        msh.faceEdges(4, ind_f2(:, end, :)) = reshape(ind_ed3(2 : end, end, :), [1, n31]);             
    end  
    if ~isPeriodicInZ    
        msh.faceEdges(1, ind_f3(:, :, end)) = reshape(ind_ed2(1 : end - 1, :, end), [1, n12]);
        msh.faceEdges(2, ind_f3(:, :, end)) = reshape(ind_ed2(2 : end, :, end), [1, n12]);
        msh.faceEdges(3, ind_f3(:, :, end)) = reshape(ind_ed1(:, 1 : end - 1, end), [1, n12]);
        msh.faceEdges(4, ind_f3(:, :, end)) = reshape(ind_ed1(:, 2 : end, end), [1, n12]);                 
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end


end














