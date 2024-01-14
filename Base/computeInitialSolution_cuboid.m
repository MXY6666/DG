% msh    : mesh of cuboid element in 3D
% u0     : cell array of exact initial solutions
% bs     : basis function set data (a struct) or just the degree of 
%          polynomial (an integer)
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function U0 = computeInitialSolution_cuboid(msh, u0, bs, layout)

if (nargin < 2)
    error('Not enough arguments')
end

if (msh.type ~= 301) && (msh.type ~= 302)
    error('Wrong mesh type')
end

if (nargin < 3) || isempty(bs)
    bs = 1;
end

if isfloat(bs) && (bs >= 1)
    k = floor(bs);
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    quad3 = GaussQuadratureRule_cube([k + 1, k + 1, k + 1], 301);
    bs    = setBasisFunctionSet_cube(quad2, quad3, 300 + k);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'cube')
        error('Wrong reference geometry for basis functions to evaluate on')
    end   
    quad3 = GaussQuadratureRule_cube(bs.neps, 301); 
    if ~(bs.elemPointsType == 301 && bs.neps(1) >= bs.deg + 1 && bs.neps(2) >= bs.deg + 1 && bs.neps(3) >= bs.deg + 1)
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        quad3 = GaussQuadratureRule_cube([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
        bs    = setBasisFunctionSet_cube(quad2, quad3, bs.type);
    end  
else
    error('Wrong argument bs')
end

if (nargin < 4) || isempty(layout)
    layout = 1;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

% Compute mass matrix at reference square
ME = computeElementMatrix_refCube([0, 0, 0], [0, 0, 0], bs);

% element center and size
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);

nv = length(u0);
U0 = zeros(bs.nb, nv * msh.nLElems);
switch layout
    case 1
        for i = 1 : nv
            U0(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems) = ME \ (bs.phitw{1} * u0{i}(ct(1, :) + 0.5 * h(1, :) .* quad3.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad3.points(:, 2), ct(3, :) + 0.5 * h(3, :) .* quad3.points(:, 3)));
        end        
    case 2
        for i = 1 : nv
            U0(:, i : nv : end) = ME \ (bs.phitw{1} * u0{i}(ct(1, :) + 0.5 * h(1, :) .* quad3.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad3.points(:, 2), ct(3, :) + 0.5 * h(3, :) .* quad3.points(:, 3)));
        end        
end


end





