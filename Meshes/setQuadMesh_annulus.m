% domain   : a 4D vector containing the radii of internal and outer arc, the
%            start and end radian, respectively
% N        : number of elements in two directions of reference domain
% type     : mesh type, 211 or 212
% maxLevel : max tree level of the mesh
% vel      : velocity vector of the flow
% ss       : sound speed
% wallType : type of the wall
% quad     : Gauss quadrature rule for a quadrilateral
% bs       : basis function set
% funy     : scaling function in radial direction
% Copyright (C) Xiangyi Meng.

function msh = setQuadMesh_annulus(domain, N, type, maxLevel, vel, ss, wallType, quad, bs, funy)

if (nargin < 2)
    error('Not enough arguments')
end

if (length(domain) ~= 4)
    error('Wrong size of argument domain')
end

if (length(N) ~= 2)
    error('Wrong size of argument N')
end

if (domain(2) <= domain(1)) || (domain(4) <= domain(3)) || ...
   ~(domain(3) >= 0 && domain(3) <= 2 * pi) || ...
   ~(domain(4) >= 0 && domain(4) <= 2 * pi)
    error('Wrong argument domain')
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

if (nargin < 10) || isempty(funy)
    funy = @(x) (exp(3 * x) - 1) / (exp(3) - 1);
end
if abs(funy(0)) > 1.0e-12 || abs(funy(1) - 1) > 1.0e-12
    error('Wrong given scaling function in radial direction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper indices
ne = prod(N);
ind_e = (1 : N(1))' + (0 : N(2) - 1) * N(1);
if (abs(domain(3)) < 1.0e-12) && (abs(domain(4) - 2 * pi) < 1.0e-12)
    ind_f2 = (1 : N(1))' + (0 : N(2)) * 2 * N(1);
    ind_f1 = [(N(1) + 1 : 2 * N(1))' + (0 : N(2) - 1) * 2 * N(1); N(1) + 1 : 2 * N(1) : 2 * ne - N(1) + 1];
    ind_n = [(1 : N(1))' + (0 : N(2)) * N(1); 1 : N(1) : ne + 1];
else
    ind_f2 = (1 : N(1))' + (0 : N(2)) * (2 * N(1) + 1);
    ind_f1 = (N(1) + 1 : 2 * N(1) + 1)' + (0 : N(2) - 1) * (2 * N(1) + 1);
    ind_n = (1 : N(1) + 1)' + (0 : N(2)) * (N(1) + 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msh          = getEmptyMesh;
msh.dm       = domain;
msh.N        = N;
msh.type     = type;
msh.maxLevel = maxLevel;
msh.bndTypes = [4 : 7, wallType];
msh.nElems   = ne;
if (abs(domain(3)) < 1.0e-12) && (abs(domain(4) - 2 * pi) < 1.0e-12)
    msh.nFaces = 2 * ne + N(1); 
    msh.nNodes = ne + N(1); 
else
    msh.nFaces = 2 * ne + N(1) + N(2);
    msh.nNodes = ne + N(1) + N(2) + 1;   
end
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems  = 1 : msh.nElems;
msh.LFaces  = 1 : msh.nFaces;
msh.nGPs    = quad.nps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the field nodes
r          = domain(1) + (domain(2) - domain(1)) * funy(linspace(0, 1, N(2) + 1));
theta      = linspace(domain(3), domain(4), N(1) + 1);
[R, Theta] = meshgrid(r, theta);
if (abs(domain(3)) < 1.0e-12) && (abs(domain(4) - 2 * pi) < 1.0e-12)
    msh.nodeCoor(1, :) = reshape(R(1 : end - 1, :) .* cos(Theta(1 : end - 1, :)), [1, msh.nNodes]);
    msh.nodeCoor(2, :) = reshape(R(1 : end - 1, :) .* sin(Theta(1 : end - 1, :)), [1, msh.nNodes]);
else
    msh.nodeCoor(1, :) = reshape(R .* cos(Theta), [1, msh.nNodes]);
    msh.nodeCoor(2, :) = reshape(R .* sin(Theta), [1, msh.nNodes]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the field elements
% (some) topology information of elements
msh.elemNodes = zeros(4, msh.nElems);
msh.elemFaces = zeros(4, msh.nElems);

msh.elemNodes(1, :) = reshape(ind_n(2 : end, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(2, :) = reshape(ind_n(1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.elemNodes(3, :) = reshape(ind_n(1 : end - 1, 2 : end), [1, msh.nElems]);
msh.elemNodes(4, :) = reshape(ind_n(2 : end, 2 : end), [1, msh.nElems]);
msh.elemFaces(1, :) = reshape(ind_f1(2 : end, :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(1 : end - 1, :), [1, msh.nElems]);
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
msh.elemIsCurved(ind_e(:, 1)) = true;

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

msh.faceNodes(1, ind_f1(1 : end - 1, :)) = reshape(ind_n(1 : end - 1, 1 : end - 1), [1, msh.nElems]);
msh.faceNodes(2, ind_f1(1 : end - 1, :)) = reshape(ind_n(1 : end - 1, 2 : end), [1, msh.nElems]);
msh.faceElems(1, ind_f1(1 : end - 1, :)) = ind_e(:)';
msh.faceElems(2, ind_f1(2 : end - 1, :)) = reshape(ind_e(1 : end - 1, :), [1, msh.nElems - N(2)]);
msh.faceNums( 1, ind_f1(1 : end - 1, :)) = 2;
msh.faceNums( 2, ind_f1(2 : end - 1, :)) = 1;

if (abs(domain(3)) < 1.0e-12) && (abs(domain(4) - 2 * pi) < 1.0e-12)
    msh.faceElems(2, ind_f1(1, :)) = ind_e(end, :);
    msh.faceNums( 2, ind_f1(1, :)) = 1;
else
    msh.faceNodes(1, ind_f1(end, :)) = ind_n(end, 2 : end);    
    msh.faceNodes(2, ind_f1(end, :)) = ind_n(end, 1 : end - 1);
    msh.faceElems(1, ind_f1(end, :)) = ind_e(end, :);
    msh.faceNums( 3, ind_f1(end, :)) = 1;     
end

msh.faceNodes(1, ind_f2) = reshape([ind_n(2 : end, 1), ind_n(1 : end - 1, 2 : end)], [1, msh.nElems + N(1)]);
msh.faceNodes(2, ind_f2) = reshape([ind_n(1 : end - 1, 1), ind_n(2 : end, 2 : end)], [1, msh.nElems + N(1)]);
msh.faceElems(1, ind_f2) = reshape([ind_e(:, 1), ind_e], [1, msh.nElems + N(1)]);
msh.faceElems(2, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 2 : end), [1, msh.nElems - N(1)]);
msh.faceNums(1, ind_f2(:, 1)) = 3;
msh.faceNums(1, ind_f2(:, 2 : end)) = 4;
msh.faceNums(2, ind_f2(:, 2 : end - 1)) = 3;

% (some) metric information of faces
faceNodes1 = msh.nodeCoor(:, msh.faceNodes(1, :));
faceNodes2 = msh.nodeCoor(:, msh.faceNodes(2, :));

msh.faceSize    = vecnorm(faceNodes2 - faceNodes1);
msh.faceNormalx = (faceNodes2(2, :) - faceNodes1(2, :)) ./ msh.faceSize;
msh.faceNormaly = (faceNodes1(1, :) - faceNodes2(1, :)) ./ msh.faceSize;
msh.faceType    = zeros(1, msh.nFaces);
if ~((abs(domain(3)) < 1.0e-12) && (abs(domain(4) - 2 * pi) < 1.0e-12))
    veln     = vel(1) * msh.faceNormalx(:, ind_f1(1, 1)) + vel(2) * msh.faceNormaly(:, ind_f1(1, 1));
    isSub    = abs(veln) < ss;
    isInflow = veln < 0;
    if isSub && isInflow
        msh.faceType(ind_f1(1, :)) = 4;
    end
    if isSub && ~isInflow
        msh.faceType(ind_f1(1, :)) = 5;
    end  
    if ~isSub && isInflow
        msh.faceType(ind_f1(1, :)) = 6;
    end    
    if ~isSub && ~isInflow
        msh.faceType(ind_f1(1, :)) = 7;
    end     
    veln     = vel(1) * msh.faceNormalx(:, ind_f1(end, 1)) + vel(2) * msh.faceNormaly(:, ind_f1(end, 1));
    isSub    = abs(veln) < ss;
    isInflow = veln < 0;    
    if isSub && isInflow
        msh.faceType(ind_f1(end, :)) = 4;
    end   
    if isSub && ~isInflow
        msh.faceType(ind_f1(end, :)) = 5;
    end  
    if ~isSub && isInflow
        msh.faceType(ind_f1(end, :)) = 6;
    end    
    if ~isSub && ~isInflow
        msh.faceType(ind_f1(end, :)) = 7;
    end    
end
msh.faceType(ind_f2(:, 1)) = wallType;
veln     = vel(1) * msh.faceNormalx(:, ind_f2(:, end)) + vel(2) * msh.faceNormaly(:, ind_f2(:, end));
isSub    = abs(veln) < ss;
isInflow = veln < 0;
msh.faceType(ind_f2(isSub  & isInflow,  end)) = 4;
msh.faceType(ind_f2(isSub  & ~isInflow, end)) = 5;
msh.faceType(ind_f2(~isSub & isInflow,  end)) = 6;
msh.faceType(ind_f2(~isSub & ~isInflow, end)) = 7;

% other useful information for computation
msh.faceJac      = 0.5 * msh.faceSize;
msh.faceIsCurved = false(1, msh.nFaces);
msh.faceIsCurved(ind_f2(:, 1)) = true;

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










