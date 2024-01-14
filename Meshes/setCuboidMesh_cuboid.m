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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% domain     : a 6D vector describing the computational domain 
% N          : 3D vector, number of elements in the x, y and z direction
% bcs        : 6D vector containing the sixth boundary types
% type = 301 : ignore the adaptive information      
% type = 302 : reserve the adaptive information
% maxLevel   : max tree level of the mesh
% funx       : scaling function in x direction
% funy       : scaling function in y direction
% funz       : scaling function in z direction
% Copyright (C) Xiangyi Meng.

function msh = setCuboidMesh_cuboid(domain, N, bcs, type, maxLevel, funx, funy, funz)
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
    type = 301;
end
if (type ~= 301) && (type ~= 302) 
    error('wrong mesh type')
end

if (nargin < 5)  || isempty(maxLevel) || (type == 301) 
    maxLevel = 0;
end
if (type == 302) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 6) || isempty(funx)
    funx = @(x)x;
end
if abs(funx(0)) > 1.0e-12 || abs(funx(1) - 1) > 1.0e-12
    error('Wrong given scaling function in x direction')
end

if (nargin < 7) || isempty(funy)
    funy = @(y)y;
end
if abs(funy(0)) > 1.0e-12 || abs(funy(1) - 1) > 1.0e-12
    error('Wrong given scaling function in y direction')
end

if (nargin < 8) || isempty(funz)
    funz = @(z)z;
end
if abs(funz(0)) > 1.0e-12 || abs(funz(1) - 1) > 1.0e-12
    error('Wrong given scaling function in z direction')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = prod(N);
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
xx = domain(1) + (domain(2) - domain(1)) * funx(linspace(0, 1, N(1) + 1));
yy = domain(3) + (domain(4) - domain(3)) * funy(linspace(0, 1, N(2) + 1));
zz = domain(5) + (domain(6) - domain(5)) * funy(linspace(0, 1, N(3) + 1));
hx = diff(xx);
hy = diff(yy); 
hz = diff(zz);
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
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems  = 1 : msh.nElems;
msh.LFaces  = 1 : msh.nFaces;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of elements
% (some) metric information of elements
msh.elemCenter = zeros(3, msh.nElems);
msh.elemLength = zeros(3, msh.nElems);
msh.elemCenter(1, :) = repmat(xx(1 : N(1)) + 0.5 * hx, [1, n23]);
msh.elemCenter(2, :) = repmat(repelem(yy(1 : N(2)) + 0.5 * hy, N(1)), [1, N(3)]);
msh.elemCenter(3, :) = repelem(zz(1 : N(3)) + 0.5 * hz, n12);
msh.elemLength(1, :) = repmat(hx, [1, n23]);
msh.elemLength(2, :) = repmat(repelem(hy, N(1)), [1, N(3)]);
msh.elemLength(3, :) = repelem(hz, n12);
msh.elemSize = prod(msh.elemLength);

% (some) topology information of elements
msh.elemFaces = zeros(6, msh.nElems);
msh.elemFaces(1, :) = reshape(ind_f1(1 : end - 1, :, :), [1, msh.nElems]);
msh.elemFaces(2, :) = reshape(ind_f1(2 : end, :, :), [1, msh.nElems]);
msh.elemFaces(3, :) = reshape(ind_f2(:, 1 : end - 1, :), [1, msh.nElems]);
msh.elemFaces(4, :) = reshape(ind_f2(:, 2 : end, :), [1, msh.nElems]);
msh.elemFaces(5, :) = reshape(ind_f3(:, :, 1 : end - 1), [1, msh.nElems]);
msh.elemFaces(6, :) = reshape(ind_f3(:, :, 2 : end), [1, msh.nElems]);

% other useful information for computation
msh.elemJac    = 0.125 * msh.elemSize;
msh.elemJxix   = 0.25 * msh.elemLength(2, :) .* msh.elemLength(3, :);
msh.elemJetay  = 0.25 * msh.elemLength(1, :) .* msh.elemLength(3, :);
msh.elemJzetaz = 0.25 * msh.elemLength(1, :) .* msh.elemLength(2, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the fields of faces
% (some) metric information of faces
msh.faceNormalx = zeros(1, msh.nFaces);
msh.faceNormaly = zeros(1, msh.nFaces);
msh.faceNormalz = zeros(1, msh.nFaces);
msh.faceSize    = zeros(1, msh.nFaces);
msh.faceType    = zeros(1, msh.nFaces);

msh.faceNormalx(ind_f1(2 : end - 1, :, :)) = 1;
msh.faceSize(   ind_f1(2 : end - 1, :, :)) = reshape(repelem(hy' * hz, N(1) - 1, 1), [1, msh.nElems - n23]);

msh.faceNormaly(ind_f2(:, 2 : end - 1, :)) = 1;
msh.faceSize(   ind_f2(:, 2 : end - 1, :)) = reshape(repelem(hx' * hz, 1, N(2) - 1), [1, msh.nElems - n31]);

msh.faceNormalz(ind_f3(:, :, 2 : end - 1)) = 1;
msh.faceSize(   ind_f3(:, :, 2 : end - 1)) = reshape(repmat(hx' * hy, [1, N(3) - 1]), [1, msh.nElems - n12]);

msh.faceNormalx(ind_f1(1, :, :)) = 1;
msh.faceSize(   ind_f1(1, :, :)) = reshape(hy' * hz, [1, n23]);
msh.faceType(   ind_f1(1, :, :)) = bcs(1);
if ~isPeriodicInX
    msh.faceNormalx(ind_f1(1, :, :))   = -1;
    msh.faceNormalx(ind_f1(end, :, :)) = 1;
    msh.faceSize(   ind_f1(end, :, :)) = reshape(hy' * hz, [1, n23]);
    msh.faceType(   ind_f1(end, :, :)) = bcs(2);
end
msh.faceNormaly(ind_f2(:, 1, :)) = 1;
msh.faceSize(   ind_f2(:, 1, :)) = reshape(hx' * hz, [1, n31]); 
msh.faceType(   ind_f2(:, 1, :)) = bcs(3);
if ~isPeriodicInY
    msh.faceNormaly(ind_f2(:, 1, :)) = -1;
    msh.faceNormaly(ind_f2(:, end, :)) = 1;
    msh.faceSize(   ind_f2(:, end, :)) = reshape(hx' * hz, [1, n31]); 
    msh.faceType(   ind_f2(:, end, :)) = bcs(4); 
end
msh.faceNormalz(ind_f3(:, :, 1)) = 1;
msh.faceSize(   ind_f3(:, :, 1)) = reshape(hx' * hy, [1, n12]);
msh.faceType(   ind_f3(:, :, 1)) = bcs(5);
if ~isPeriodicInZ
    msh.faceNormalz(ind_f3(:, :, 1)) = -1;
    msh.faceNormalz(ind_f3(:, :, end)) = 1;
    msh.faceSize(   ind_f3(:, :, end)) = reshape(hx' * hy, [1, n12]);
    msh.faceType(   ind_f3(:, :, end)) = bcs(6);    
end

% (some) topology information of faces
msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums  = zeros(2, msh.nFaces); 

msh.faceElems(1, ind_f1(2 : end - 1, :, :)) = reshape(ind_e(1 : end - 1, :, :), [1, msh.nElems - n23]);
msh.faceElems(2, ind_f1(2 : end - 1, :, :)) = reshape(ind_e(2 : end, :, :), [1, msh.nElems - n23]);
msh.faceNums( 1, ind_f1(2 : end - 1, :, :)) = 2;
msh.faceNums( 2, ind_f1(2 : end - 1, :, :)) = 1;

msh.faceElems(1, ind_f2(:, 2 : end - 1, :)) = reshape(ind_e(:, 1 : end - 1, :), [1, msh.nElems - n31]);
msh.faceElems(2, ind_f2(:, 2 : end - 1, :)) = reshape(ind_e(:, 2 : end, :), [1, msh.nElems - n31]);
msh.faceNums( 1, ind_f2(:, 2 : end - 1, :)) = 4;
msh.faceNums( 2, ind_f2(:, 2 : end - 1, :)) = 3;

msh.faceElems(1, ind_f3(:, :, 2 : end - 1)) = reshape(ind_e(:, :, 1 : end - 1), [1, msh.nElems - n12]);
msh.faceElems(2, ind_f3(:, :, 2 : end - 1)) = reshape(ind_e(:, :, 2 : end), [1, msh.nElems - n12]);
msh.faceNums( 1, ind_f3(:, :, 2 : end - 1)) = 6;
msh.faceNums( 2, ind_f3(:, :, 2 : end - 1)) = 5;

if isPeriodicInX 
    msh.faceElems(1, ind_f1(1, :, :)) = reshape(ind_e(end, :, :), [1, n23]);
    msh.faceElems(2, ind_f1(1, :, :)) = reshape(ind_e(1, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(1, :, :)) = 2;
    msh.faceNums( 2, ind_f1(1, :, :)) = 1;    
else
    msh.faceElems(1, ind_f1(1, :, :))   = reshape(ind_e(1, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(1, :, :))   = 1;  
    msh.faceElems(1, ind_f1(end, :, :)) = reshape(ind_e(end, :, :), [1, n23]);
    msh.faceNums( 1, ind_f1(end, :, :)) = 2;    
end
if isPeriodicInY
    msh.faceElems(1, ind_f2(:, 1, :)) = reshape(ind_e(:, end, :), [1, n31]);
    msh.faceElems(2, ind_f2(:, 1, :)) = reshape(ind_e(:, 1, :), [1, n31]);
    msh.faceNums( 1, ind_f2(:, 1, :)) = 4;
    msh.faceNums( 2, ind_f2(:, 1, :)) = 3;    
else
    msh.faceElems(1, ind_f2(:, 1, :))   = reshape(ind_e(:, 1, :), [1, n31]);
    msh.faceNums( 1, ind_f2(:, 1, :))   = 3;
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
    msh.faceElems(1, ind_f3(:, :, end)) = reshape(ind_e(:, :, end), [1, n12]);
    msh.faceNums( 1, ind_f3(:, :, end)) = 6;     
end

% other useful information for computation
msh.faceJac = 0.25 * msh.faceSize;

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
if (type == 302)
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
end

end

















