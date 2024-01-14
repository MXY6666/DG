% The node number is indexed like this
%            4---------3
%            |         |
%            |         |
%            |         |
%            |         |
%            |         |
%            1---------2
% The face number is indexed like this
%            o----4----o
%            |         |
%            |         |
%            1         2
%            |         |
%            |         |
%            o----3----o
%
%            o--12--o--11--o
%            |             |
%            05            08
%            |             |
%            o             o
%            |             |
%            06            07
%            |             |
%            o--09--o--10--o
% 
% See getEmptyMesh.m for detailed fields of mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain     : a 4D vector describing the computational domain 
% N          : 2D vector, number of elements in the x and y direction
% bcs        : 1 * 4 or 4 * 1 vector containing the four boundary types
% type = 211 : ignore the adaptive information   
% type = 212 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% pert       : a 2D vector containing amplitude of nodes perturbation in x and y
%              diretion, respectively
% quad       : Gauss quadrature rule for a quadrilateral
% bs         : basis function set
% Copyright (C) Xiangyi Meng.

function msh = setQuadMesh_rect(domain, N, bcs, type, maxLevel, pert, quad, bs)

if (nargin < 2)
    error('Not enough arguments')
end

if (length(domain) ~= 4)
    error('Wrong size of argument domain')
end

if (length(N) ~= 2)
    error('Wrong size of argument N')
end

if (domain(2) <= domain(1)) || (domain(4) <= domain(3))
    error('Wrong argument domain')
end

if (nargin < 3)  || isempty(bcs)
    bcs = ones(1, 4);
end
if (length(bcs) ~= 4)
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

if (nargin < 4)  || isempty(type)
    type = 211;
end
if (type ~= 211) && (type ~= 212) 
    error('wrong mesh type')
end

if (nargin < 5)  || isempty(maxLevel) || (type == 211) 
    maxLevel = 0;
end
if (type == 212) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 6)  || isempty(pert)
    pert = [0.3, 0.3];
end
if (length(pert) ~= 2 || abs(pert(1)) >= 0.5 || abs(pert(2)) >= 0.5)
    error('Wrong amplitude of nodes perturbation')
end

if (nargin < 7) || isempty(quad)
    quad = GaussQuadratureRule_quad([2, 2], 201);
end

if (nargin < 8) || isempty(bs)
    quad1 = GaussQuadratureRule_line(2, 101);
    bs = setBasisFunctionSet_square(quad1, quad, 201, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = prod(N);
nn = ne + sum(N) + 1;
ind_e = (1 : N(1))' + (0 : N(2) - 1) * N(1);
if isPeriodicInX
    if isPeriodicInY
        ind_f2 = [(1 : N(1))' + (0 : N(2) - 1) * 2 * N(1), (1 : N(1))'];
        ind_f1 = [(N(1) + 1 : 2 * N(1))' + (0 : N(2) - 1) * 2 * N(1); N(1) + 1 : 2 * N(1) : 2 * ne - N(1) + 1];
    else
        ind_f2 = (1 : N(1))' + (0 : N(2)) * 2 * N(1);
        ind_f1 = [(N(1) + 1 : 2 * N(1))' + (0 : N(2) - 1) * 2 * N(1); N(1) + 1 : 2 * N(1) : 2 * ne - N(1) + 1];
    end
else
    if isPeriodicInY
        ind_f2 = [(1 : N(1))' + (0 : N(2) - 1) * (2 * N(1) + 1), (1 : N(1))'];
        ind_f1 = (N(1) + 1 : 2 * N(1) + 1)' + (0 : N(2) - 1) * (2 * N(1) + 1);
    else
        ind_f2 = (1 : N(1))' + (0 : N(2)) * (2 * N(1) + 1);
        ind_f1 = (N(1) + 1 : 2 * N(1) + 1)' + (0 : N(2) - 1) * (2 * N(1) + 1);
    end
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
msh.nFaces   = 2 * ne + sum(N); 
if isPeriodicInX
    msh.nFaces = msh.nFaces - N(2);
end
if isPeriodicInY
    msh.nFaces = msh.nFaces - N(1);
end
msh.nNodes  = nn;
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems  = 1 : msh.nElems;
msh.LFaces  = 1 : msh.nFaces;
msh.nGPs    = quad.nps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of nodes
msh.nodeCoor = zeros(2, msh.nNodes);
msh.nodeCoor(1, :) = repmat(linspace(domain(1), domain(2), N(1) + 1), [1, N(2) + 1]);
msh.nodeCoor(2, :) = repelem(linspace(domain(3), domain(4), N(2) + 1), N(1) + 1);

% Perturb the nodes
if (abs(pert(1)) > 0)
    rands = pert(1) * (2 * rand(N(1) - 1, N(2) + 1) - 1);
    if isPeriodicInY
        rands(:, end) = rands(:, 1);
    end
    hx = (domain(2) - domain(1)) / N(1);
    msh.nodeCoor(1, ind_n(2 : end - 1, :)) = msh.nodeCoor(1, ind_n(2 : end - 1, :)) + hx * rands(:)';
end
if (abs(pert(2)) > 0)
    rands = pert(2) * (2 * rand(N(1) + 1, N(2) - 1) - 1);
    if isPeriodicInX
        rands(end, :) = rands(1, :);
    end
    hy = (domain(4) - domain(3)) / N(2);
    msh.nodeCoor(2, ind_n(:, 2 : end - 1)) = msh.nodeCoor(2, ind_n(:, 2 : end - 1)) + hy * rands(:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of elements
% (some) topology information of elements
msh.elemNodes = zeros(4, msh.nElems);
msh.elemFaces = zeros(4, msh.nElems);

msh.elemNodes(1, :) = reshape(ind_n(1 : N(1), 1 : N(2)), [1, msh.nElems]);
msh.elemNodes(2, :) = reshape(ind_n(2 : end, 1 : N(2)), [1, msh.nElems]);
msh.elemNodes(3, :) = reshape(ind_n(2 : end, 2 : end), [1, msh.nElems]);
msh.elemNodes(4, :) = reshape(ind_n(1 : N(1), 2 : end), [1, msh.nElems]);
msh.elemFaces(1, :) = reshape(ind_f1(1 : N(1), :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(2 : end, :), [1, msh.nElems]);
msh.elemFaces(3, :) = reshape(ind_f2(:, 1 : N(2)), [1, msh.nElems]);
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

% Store inverse of mass matrix in each element
massMat        = bs.phitw{1} * (repmat(bs.phi{1}, [1, msh.nElems]) .* repelem(msh.elemJac, 1, bs.nb));
massMatInv     = cellfun(@inv, mat2cell(massMat, bs.nb, bs.nb * ones(1, msh.nElems)), 'UniformOutput', false);
msh.massMatInv = cellfun(@mychop, massMatInv, num2cell(1.0e-12 * ones(1, msh.nElems)), 'UniformOutput', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of faces
% (some) topology information of faces
msh.faceNodes = zeros(2, msh.nFaces);
msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums  = zeros(2, msh.nFaces); 

msh.faceNodes(1, ind_f1(2 : end - 1, :)) = reshape(ind_n(2 : end - 1, 1 : end - 1), [1, msh.nElems - N(2)]);
msh.faceNodes(2, ind_f1(2 : end - 1, :)) = reshape(ind_n(2 : end - 1, 2 : end), [1, msh.nElems - N(2)]);
msh.faceElems(1, ind_f1(2 : end - 1, :)) = reshape(ind_e(1 : end - 1, :), [1, msh.nElems - N(2)]);
msh.faceElems(2, ind_f1(2 : end - 1, :)) = reshape(ind_e(2 : end, :), [1, msh.nElems - N(2)]);
msh.faceNums( 1, ind_f1(2 : end - 1, :)) = 2;
msh.faceNums( 2, ind_f1(2 : end - 1, :)) = 1;

msh.faceNodes(1, ind_f2(:, 2 : end - 1)) = reshape(ind_n(2 : end, 2 : end - 1), [1, msh.nElems - N(1)]);
msh.faceNodes(2, ind_f2(:, 2 : end - 1)) = reshape(ind_n(1 : end - 1, 2 : end - 1), [1, msh.nElems - N(1)]);
msh.faceElems(1, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 1 : end - 1), [1, msh.nElems - N(1)]);
msh.faceElems(2, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 2 : end), [1, msh.nElems - N(1)]);
msh.faceNums( 1, ind_f2(:, 2 : end - 1)) = 4;
msh.faceNums( 2, ind_f2(:, 2 : end - 1)) = 3;

if isPeriodicInX 
    msh.faceNodes(1, ind_f1(1, :)) = ind_n(1, 1 : end - 1);
    msh.faceNodes(2, ind_f1(1, :)) = ind_n(1, 2 : end);
    msh.faceElems(1, ind_f1(1, :)) = ind_e(end, :);
    msh.faceElems(2, ind_f1(1, :)) = ind_e(1, :);
    msh.faceNums( 1, ind_f1(1, :)) = 2;
    msh.faceNums( 2, ind_f1(1, :)) = 1;    
else
    msh.faceNodes(1, ind_f1(1, :))   = ind_n(1, 2 : end);
    msh.faceNodes(2, ind_f1(1, :))   = ind_n(1, 1 : end - 1);    
    msh.faceElems(1, ind_f1(1, :))   = ind_e(1, :);
    msh.faceNums( 1, ind_f1(1, :))   = 1;
    msh.faceNodes(1, ind_f1(end, :)) = ind_n(end, 1 : end - 1);
    msh.faceNodes(2, ind_f1(end, :)) = ind_n(end, 2 : end);    
    msh.faceElems(1, ind_f1(end, :)) = ind_e(end, :);
    msh.faceNums( 1, ind_f1(end, :)) = 2;    
end
if isPeriodicInY
    msh.faceNodes(1, ind_f2(:, 1)) = ind_n(2 : end, 1);
    msh.faceNodes(2, ind_f2(:, 1)) = ind_n(1 : end - 1, 1);    
    msh.faceElems(1, ind_f2(:, 1)) = ind_e(:, end)';
    msh.faceElems(2, ind_f2(:, 1)) = ind_e(:, 1)';
    msh.faceNums( 1, ind_f2(:, 1)) = 4;
    msh.faceNums( 2, ind_f2(:, 1)) = 3;    
else
    msh.faceNodes(1, ind_f2(:, 1))   = ind_n(1 : end - 1, 1);
    msh.faceNodes(2, ind_f2(:, 1))   = ind_n(2 : end, 1);    
    msh.faceElems(1, ind_f2(:, 1))   = ind_e(:, 1)';
    msh.faceNums( 1, ind_f2(:, 1))   = 3;
    msh.faceNodes(1, ind_f2(:, end)) = ind_n(2 : end, end);
    msh.faceNodes(2, ind_f2(:, end)) = ind_n(1 : end - 1, end);    
    msh.faceElems(1, ind_f2(:, end)) = ind_e(:, end)';
    msh.faceNums( 1, ind_f2(:, end)) = 4;   
end

% (some) metric information of faces
faceNodes1 = msh.nodeCoor(:, msh.faceNodes(1, :));
faceNodes2 = msh.nodeCoor(:, msh.faceNodes(2, :));

msh.faceSize = vecnorm(faceNodes2 - faceNodes1);
msh.faceNormalx = (faceNodes2(2, :) - faceNodes1(2, :)) ./ msh.faceSize;
msh.faceNormaly = (faceNodes1(1, :) - faceNodes2(1, :)) ./ msh.faceSize;
msh.faceType = zeros(1, msh.nFaces);
msh.faceType(ind_f1(1, :)) = bcs(1);
msh.faceType(ind_f2(:, 1)) = bcs(3);
if ~isPeriodicInX 
    msh.faceType(ind_f1(end, :)) = bcs(2);
end
if ~isPeriodicInY
    msh.faceType(ind_f2(:, end)) = bcs(4); 
end

% other useful information for computation
msh.faceJac      = 0.5 * msh.faceSize;
msh.faceIsCurved = false(1, msh.nFaces);

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

