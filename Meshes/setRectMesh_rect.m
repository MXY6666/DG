% The face number is indexed like this
%            o----4----o
%            |         |
%            |         |
%            1         2
%            |         |
%            |         |
%            o----3----o
%
%            o--11--o--12--o
%            |             |
%            06            08
%            |             |
%            o             o
%            |             |
%            05            07
%            |             |
%            o--09--o--10--o
%
% See getEmptyMesh.m for detailed fields of mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain     : a 4D vector describing the computational domain 
% N          : 2D vector, number of elements in the x and y direction
% bcs        : 1 * 4 or 4 * 1 vector containing the four boundary types
% type = 201 : ignore the adaptive information         
% type = 202 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% funx       : scaling function in x direction
% funy       : scaling function in y direction
% Copyright (C) Xiangyi Meng.

function msh = setRectMesh_rect(domain, N, bcs, type, maxLevel, funx, funy)

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
    type = 201;
end
if (type ~= 201) && (type ~= 202) 
    error('wrong mesh type')
end

if (nargin < 5)  || isempty(maxLevel) || (type == 201) 
    maxLevel = 0;
end
if (type == 202) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 6) || isempty(funx)
    funx = @(x)x;
    % funx = @(x)(exp(3 * x) - 1) / (exp(3) - 1);
    % funx = @(x)0.5 * (1 + tanh(3 * (2 * x - 1)) / tanh(3));
end
if abs(funx(0)) > 1.0e-12 || abs(funx(1) - 1) > 1.0e-12
    error('Wrong given scaling function in x direction')
end

if (nargin < 7) || isempty(funy)
    funy = @(y)y;
    % funy = @(y)(exp(3 * y) - 1) / (exp(3) - 1);
    % funy = @(y)0.5 * (1 + tanh(3 * (2 * y - 1)) / tanh(3));
end
if abs(funy(0)) > 1.0e-12 || abs(funy(1) - 1) > 1.0e-12
    error('Wrong given scaling function in y direction')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = prod(N);
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
xx = domain(1) + (domain(2) - domain(1)) * funx(linspace(0, 1, N(1) + 1));
yy = domain(3) + (domain(4) - domain(3)) * funy(linspace(0, 1, N(2) + 1));
hx = diff(xx);
hy = diff(yy);  
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
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems  = 1 : msh.nElems;
msh.LFaces  = 1 : msh.nFaces;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of elements
% (some) metric information of elements
msh.elemCenter = [repmat(xx(1 : N(1)) + 0.5 * hx, [1, N(2)]); repelem(yy(1 : N(2)) + 0.5 * hy, N(1))];
msh.elemLength = [repmat(hx, [1, N(2)]); repelem(hy, N(1))];
msh.elemSize = prod(msh.elemLength);

% (some) topology information of elements
msh.elemFaces = zeros(4, msh.nElems);
msh.elemFaces(1, :) = reshape(ind_f1(1 : N(1), :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(2 : end, :), [1, msh.nElems]);
msh.elemFaces(3, :) = reshape(ind_f2(:, 1 : N(2)), [1, msh.nElems]);
msh.elemFaces(4, :) = reshape(ind_f2(:, 2 : end), [1, msh.nElems]);

% other useful information for computation
msh.elemJac   = 0.25 * msh.elemSize;
msh.elemJxix  = 0.5 * msh.elemLength(2, :);
msh.elemJetay = 0.5 * msh.elemLength(1, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of faces
% (some) metric information of faces
msh.faceNormalx = zeros(1, msh.nFaces);
msh.faceNormaly = zeros(1, msh.nFaces);
msh.faceSize    = zeros(1, msh.nFaces);
msh.faceType    = zeros(1, msh.nFaces);

msh.faceNormalx(ind_f1(2 : end - 1, :)) = 1;
msh.faceSize(   ind_f1(2 : end - 1, :)) = repelem(hy, N(1) - 1);

msh.faceNormaly(ind_f2(:, 2 : end - 1)) = 1;
msh.faceSize(   ind_f2(:, 2 : end - 1)) = repmat(hx, [1, N(2) - 1]);

msh.faceNormalx(ind_f1(1, :)) = 1;
msh.faceSize(   ind_f1(1, :)) = hy;
msh.faceType(   ind_f1(1, :)) = bcs(1);
if ~isPeriodicInX
    msh.faceNormalx(ind_f1(1, :))   = -1;
    msh.faceNormalx(ind_f1(end, :)) = 1;
    msh.faceSize(   ind_f1(end, :)) = hy;  
    msh.faceType(   ind_f1(end, :)) = bcs(2);
end
msh.faceNormaly(ind_f2(:, 1)) = 1;
msh.faceSize(   ind_f2(:, 1)) = hx; 
msh.faceType(   ind_f2(:, 1)) = bcs(3);
if ~isPeriodicInY
    msh.faceNormaly(ind_f2(:, 1)) = -1;
    msh.faceNormaly(ind_f2(:, end)) = 1;
    msh.faceSize(   ind_f2(:, end)) = hx; 
    msh.faceType(   ind_f2(:, end)) = bcs(4); 
end

% (some) topology information of faces
msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums  = zeros(2, msh.nFaces); 

msh.faceElems(1, ind_f1(2 : end - 1, :)) = reshape(ind_e(1 : end - 1, :), [1, msh.nElems - N(2)]);
msh.faceElems(2, ind_f1(2 : end - 1, :)) = reshape(ind_e(2 : end, :), [1, msh.nElems - N(2)]);
msh.faceNums( 1, ind_f1(2 : end - 1, :)) = 2;
msh.faceNums( 2, ind_f1(2 : end - 1, :)) = 1;

msh.faceElems(1, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 1 : end - 1), [1, msh.nElems - N(1)]);
msh.faceElems(2, ind_f2(:, 2 : end - 1)) = reshape(ind_e(:, 2 : end), [1, msh.nElems - N(1)]);
msh.faceNums( 1, ind_f2(:, 2 : end - 1)) = 4;
msh.faceNums( 2, ind_f2(:, 2 : end - 1)) = 3;

if isPeriodicInX 
    msh.faceElems(1, ind_f1(1, :)) = ind_e(end, :);
    msh.faceElems(2, ind_f1(1, :)) = ind_e(1, :);
    msh.faceNums( 1, ind_f1(1, :)) = 2;
    msh.faceNums( 2, ind_f1(1, :)) = 1;    
else
    msh.faceElems(1, ind_f1(1, :))   = ind_e(1, :);
    msh.faceNums( 1, ind_f1(1, :))   = 1;  
    msh.faceElems(1, ind_f1(end, :)) = ind_e(end, :);
    msh.faceNums( 1, ind_f1(end, :)) = 2;    
end
if isPeriodicInY
    msh.faceElems(1, ind_f2(:, 1)) = ind_e(:, end)';
    msh.faceElems(2, ind_f2(:, 1)) = ind_e(:, 1)';
    msh.faceNums( 1, ind_f2(:, 1)) = 4;
    msh.faceNums( 2, ind_f2(:, 1)) = 3;    
else
    msh.faceElems(1, ind_f2(:, 1))   = ind_e(:, 1)';
    msh.faceNums( 1, ind_f2(:, 1))   = 3;
    msh.faceElems(1, ind_f2(:, end)) = ind_e(:, end)';
    msh.faceNums( 1, ind_f2(:, end)) = 4;    
end

% other useful information for computation
msh.faceJac = 0.5 * msh.faceSize;

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
if (type == 202)
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


