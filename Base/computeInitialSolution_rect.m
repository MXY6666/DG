% msh    : mesh of rectangular element in 2D
% u0     : cell array of exact initial solutions
% bs     : basis function set data (a struct) or just the degree of 
%          polynomial(an integer)
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function U0 = computeInitialSolution_rect(msh, u0, bs, layout)

if (nargin < 2)
    error('Not enough arguments')
end

if (msh.type ~= 201) && (msh.type ~= 202)
    error('Wrong mesh type')
end

if (nargin < 3) || isempty(bs)
    bs = 1;
end

if isfloat(bs) && (bs >= 1)
    k     = floor(bs);
    quad1 = GaussQuadratureRule_line(k + 1, 101); 
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    bs    = setBasisFunctionSet_square(quad1, quad2, 200 + k);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'square')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    quad2 = GaussQuadratureRule_square(bs.neps, bs.elemPointsType);
    if (bs.elemPointsType == 201) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 1)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
    end
    if (bs.elemPointsType == 202) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 2)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 2, bs.deg + 2], 202);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
    end       
    if (bs.elemPointsType == 203) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 2)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 2], 203);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
    end     
    if (bs.elemPointsType == 204) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 1)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 2, bs.deg + 1], 204);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
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
ME = computeElementMatrix_refSquare([0, 0], [0, 0], bs);

% element center and size
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);

nv = length(u0);
U0 = zeros(bs.nb, nv * msh.nLElems);
switch layout
    case 1
        for i = 1 : nv
            U0(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems) = ME \ (bs.phitw{1} * u0{i}(ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2)));
        end        
    case 2
        for i = 1 : nv
            U0(:, i : nv : end) = ME \ (bs.phitw{1} * u0{i}(ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2)));
        end        
end

end
