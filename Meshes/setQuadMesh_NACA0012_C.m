% domain   : a 4D vector containing the coordinates of the leading edge, 
%            radius of the left semi-circle which centered at the leading 
%            edge of the airfoil, and the length of right rectangle
% N        : a 3D vector containing number of elements in two directions of 
%            reference domain and number of elements near the airfoil
% type     : mesh type, 211 or 212
% maxLevel : max tree level of the mesh
% vel      : velocity vector of the flow
% ss       : sound speed
% wallType : type of the wall
% quad     : Gauss quadrature rule for a quadrilateral
% bs       : basis function set
% refine   : true or false, indicating if refine the elements near the solid 
%            wall and ones behind the trailing edge
% n        : if the choice 'refine' still cannot satisfy your needs, you can
%           provide n which means the ith row of elements will be averagely 
%           split into n(i) rows
% Copyright (C) Xiangyi Meng.

function msh = setQuadMesh_NACA0012_C(domain, N, type, maxLevel, vel, ss, wallType, quad, bs, refine, n)

if (nargin < 2)
    error('Not enough arguments')
end

if (length(domain) ~= 4) || (domain(4) <= 1) 
    error('Wrong argument domain')
end

if (length(N) ~= 3) || (mod(N(1) - N(3), 2) ~= 0)
    error('Wrong argument N')
end

if (nargin < 3)  || isempty(type)
    type = 211;
end
if (type ~= 211) && (type ~= 212) 
    error('wrong mesh type')
end

if (nargin < 4)  || isempty(maxLevel) || (type == 211) 
    maxLevel = 0;
end
if (type == 212) && (maxLevel == 0)
    maxLevel = 3;  
end

if (nargin < 5)  || isempty(vel)
    vel = [1; 0];
end

if (nargin < 6)  || isempty(ss)
    ss = 1;
end

if (nargin < 7)  || isempty(wallType)
    wallType = 8;
end

if (nargin < 8) || isempty(quad)
    quad = GaussQuadratureRule_quad([2, 2], 201);
end

if (nargin < 9) || isempty(bs)
    quad1 = GaussQuadratureRule_line(2, 101);
    bs = setBasisFunctionSet_square(quad1, quad, 201, 2);
end

if (nargin < 10)  || isempty(refine) 
    refine = false;
end

if (nargin < 11)
    n = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper indices
ln = length(n); sn = sum(n);
N(2) = N(2) - ln + sn;
cut1 = (N(1) - N(3)) / 2; cut2 = (N(1) + N(3)) / 2;
ne = prod(N(1 : 2));
ind_e = reshape(1 : ne, N(1 : 2));
ind_f2 = (1 : N(1))' + (0 : N(2)) * (2 * N(1) + 1);
ind_f1 = (N(1) + 1 : 2 * N(1) + 1)' + (0 : N(2) - 1) * (2 * N(1) + 1);
ind_n = reshape(1 : (N(1) + 1) * (N(2) + 1), N(1 : 2) + 1);
ind_f2(cut2 + 1 : end, 1) = ind_f2(cut1 : -1 : 1, 1);
ind_f2(:, 2 : end) = ind_f2(:, 2 : end) - cut1;
ind_f1 = ind_f1 - cut1;
ind_n(cut2 + 1 : end, 1) = ind_n(cut1 + 1 : -1 : 1, 1);
ind_n(:, 2 : end) = ind_n(:, 2 : end) - (cut1 + 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msh          = getEmptyMesh;
msh.dm       = domain;
msh.N        = N;
msh.type     = type;
msh.maxLevel = maxLevel;
msh.bndTypes = [4 : 7, wallType];
msh.nElems   = ne;
msh.nFaces   = 2 * ne + N(2) + cut2;
msh.nNodes   = ne + N(2) + cut2;
msh.nLElems  = msh.nElems;
msh.nLFaces  = msh.nFaces;
msh.LElems   = 1 : msh.nElems;
msh.LFaces   = 1 : msh.nFaces;
msh.nGPs     = quad.nps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the field nodes
a = zeros(2, 4); 
a(1, 1) = domain(3); a(1, 2) = domain(4); a(2, 1) = N(3) / N(1);
if refine
    it = [0, 0, 2, 2];
else
    it = zeros(1, 4);
end
[x, y] = tensor(N(1), N(2) - sn + ln, 21, a, it, 14);

xnew = zeros(N(1) + 1, N(2) + 1);
ynew = zeros(N(1) + 1, N(2) + 1);
xnew(:, sn + 1 : end) = x(:, ln + 1 : end);
ynew(:, sn + 1 : end) = y(:, ln + 1 : end);
s = 1;
for j = 1 : ln
    dx = (x(:, j + 1) - x(:, j)) / n(j);
    dy = (y(:, j + 1) - y(:, j)) / n(j);
    xnew(:, s : s + n(j) - 1) = x(:, j) + (0 : n(j) - 1) .* dx;
    ynew(:, s : s + n(j) - 1) = y(:, j) + (0 : n(j) - 1) .* dy;
    s = s + n(j);
end
msh.nodeCoor(1, :) = [xnew(1 : cut2), xnew(N(1) + 2 : end)] + domain(1); 
msh.nodeCoor(2, :) = [ynew(1 : cut2), ynew(N(1) + 2 : end)] + domain(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the field elements
% (some) topology information of elements
msh.elemNodes = zeros(4, msh.nElems);
msh.elemFaces = zeros(4, msh.nElems);

msh.elemNodes(1, :) = reshape(ind_n(1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(2, :) = reshape(ind_n(2 : end, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(3, :) = reshape(ind_n(2 : end, 2 : end), [1, msh.nElems]);
msh.elemNodes(4, :) = reshape(ind_n(1 : end - 1, 2 : end), [1, msh.nElems]);
msh.elemFaces(1, :) = reshape(ind_f1(1 : end - 1, :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(2 : end, :), [1, msh.nElems]);
msh.elemFaces(3, :) = reshape(ind_f2(:, 1 : end - 1), [1, msh.nElems]);
msh.elemFaces(4, :) = reshape(ind_f2(:, 2 : end), [1, msh.nElems]);
   
% (some) metric information of elements
x = reshape(msh.nodeCoor(1, msh.elemNodes), [4, msh.nElems]);
y = reshape(msh.nodeCoor(2, msh.elemNodes), [4, msh.nElems]);
msh.elemSize = 0.5 * ((x(1, :) - x(3, :)) .* (y(2, :) - y(4, :)) - (x(2, :) - x(4, :)) .* (y(1, :) - y(3, :)));
msh.elemDiam = sqrt(msh.elemSize);

% other useful information for computation
xxi  = quad.N_xi  * x;
xeta = quad.N_eta * x;
yxi  = quad.N_xi  * y;
yeta = quad.N_eta * y;

msh.elemGPx = quad.N * x;
msh.elemGPy = quad.N * y;
msh.elemJac = xxi .* yeta - xeta .* yxi;
if any(msh.elemJac < 0)
    error('Negative Jacobian')
end
msh.elemJxix  = yeta; 
msh.elemJxiy  = -xeta; 
msh.elemJetax = -yxi; 
msh.elemJetay = xxi; 
msh.elemIsCurved = false(1, msh.nElems);
msh.elemIsCurved(ind_e(cut1 + 1 : cut2, 1)) = true;

% Store inverse of mass matrix in each element
massMat        = bs.phitw{1} * (repmat(bs.phi{1}, [1, msh.nElems]) .* repelem(msh.elemJac, 1, bs.nb));
massMatInv     = cellfun(@inv, mat2cell(massMat, bs.nb, bs.nb * ones(1, msh.nElems)), 'UniformOutput', false);
msh.massMatInv = cellfun(@mychop, massMatInv, num2cell(1.0e-12 * ones(1, msh.nElems)), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the field faces
% (some) topology information of faces
msh.faceNodes = zeros(2, msh.nFaces);
msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums  = zeros(2, msh.nFaces); 

% Initialize the field faces
msh.faceNodes(1, 1 : cut2) = 1 : cut2;
msh.faceNodes(1, ind_f2(:, 2 : end)) = reshape(ind_n(2 : end, 2 : end), [1, msh.nElems]);
msh.faceNodes(2, 1 : cut2) = ind_n(2 : cut2 + 1, 1)';
msh.faceNodes(2, ind_f2(:, 2 : end)) = reshape(ind_n(1 : end - 1, 2 : end), [1, msh.nElems]);
msh.faceElems(1, 1 : cut2) = 1 : cut2;
msh.faceElems(1, ind_f2(:, 2 : end)) = ind_e(:)';
msh.faceElems(2, 1 : cut1) = N(1) : -1 : cut2 + 1;
msh.faceElems(2, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 2 : end), [1, msh.nElems - N(1)]);
msh.faceNums( 1, 1 : cut2) = 3;
msh.faceNums( 1, ind_f2(:, 2 : end)) = 4;
msh.faceNums( 2, 1 : cut1) = 3;
msh.faceNums( 2, ind_f2(:, 2 : end - 1)) = 3;

msh.faceNodes(1, ind_f1) = reshape([ind_n(1, 2 : end); ind_n(2 : end, 1 : end - 1)], [1, msh.nElems + N(2)]);
msh.faceNodes(2, ind_f1) = reshape([ind_n(1, 1 : end - 1); ind_n(2 : end, 2 : end)], [1, msh.nElems + N(2)]);
msh.faceElems(1, ind_f1) = reshape([ind_e(1, :); ind_e], [1, msh.nElems + N(2)]);
msh.faceElems(2, ind_f1(2 : end - 1, :)) = reshape(ind_e(2 : end, :), [1, msh.nElems - N(2)]);
msh.faceNums( 1, ind_f1(1, :)) = 1;
msh.faceNums( 1, ind_f1(2 : end, :)) = 2;
msh.faceNums( 2, ind_f1(2 : end - 1, :)) = 1;

% (some) metric information of faces
faceNodes1 = msh.nodeCoor(:, msh.faceNodes(1, :));
faceNodes2 = msh.nodeCoor(:, msh.faceNodes(2, :));

msh.faceSize    = vecnorm(faceNodes2 - faceNodes1);
msh.faceNormalx = (faceNodes2(2, :) - faceNodes1(2, :)) ./ msh.faceSize;
msh.faceNormaly = (faceNodes1(1, :) - faceNodes2(1, :)) ./ msh.faceSize;
msh.faceType    = zeros(1, msh.nFaces);

msh.faceType(cut1 + 1 : cut2) = wallType;
veln     = vel(1) * msh.faceNormalx(:, ind_f2(:, end)) + vel(2) * msh.faceNormaly(:, ind_f2(:, end));
isSub    = abs(veln) < ss;
isInflow = veln < 0;
msh.faceType(ind_f2(isSub  & isInflow,  end)) = 4;
msh.faceType(ind_f2(isSub  & ~isInflow, end)) = 5;
msh.faceType(ind_f2(~isSub & isInflow,  end)) = 6;
msh.faceType(ind_f2(~isSub & ~isInflow, end)) = 7;
veln  = vel(1) * msh.faceNormalx(:, ind_f1(1, 1)) + vel(2) * msh.faceNormaly(:, ind_f1(1, 1));
isSub = abs(veln) < ss;
if isSub
    msh.faceType(ind_f1([1, end], :)) = 5;
else
    msh.faceType(ind_f1([1, end], :)) = 7;
end

% other useful information for computation
msh.faceJac      = 0.5 * msh.faceSize;
msh.faceIsCurved = false(1, msh.nFaces);
msh.faceIsCurved(ind_f2(cut1 + 1 : cut2, 1)) = true;

% classify internal faces and boundary faces
msh.intLFaces  = find(msh.faceType == 0);
msh.nIntLFaces = length(msh.intLFaces);

msh.bndLFaces  = cell(1, length(msh.bndTypes));
msh.nBndLFaces = zeros(1, length(msh.bndTypes));
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i}  = find(msh.faceType == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (type == 212)
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

end













