% msh    : mesh of line element in 1D
% u0     : cell array of exact initial solutions
% bs     : basis function set data (a struct) or just the degree of 
%          polynomial (an integer)
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function U0 = computeInitialSolution_line(msh, u0, bs, layout)

if (nargin < 2)
    error('Not enough arguments')
end

if (msh.type ~= 101) && (msh.type ~= 102)
    error('Wrong mesh type')
end

if (nargin < 3) || isempty(bs)
    bs = 1;
end

if isfloat(bs) && (bs >= 1)
    k    = floor(bs); 
    quad = GaussQuadratureRule_line(k + 1, 101); 
    bs   = setBasisFunctionSet_line(quad, 200 + k);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'line')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    quad = GaussQuadratureRule_line(bs.nep, bs.elemPointsType);
    if (bs.elemPointsType == 101) && (bs.nep < bs.deg + 1)                                                     
        quad = GaussQuadratureRule_line(bs.deg + 1, 101);
        bs   = setBasisFunctionSet_line(quad, bs.type);
    end
    if (bs.elemPointsType == 102) && (bs.nep < bs.deg + 2)
        quad = GaussQuadratureRule_line(bs.deg + 2, 102);
        bs   = setBasisFunctionSet_line(quad, bs.type);                                                                            
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

% Compute mass matrix at reference line
ME = computeElementMatrix_refLine(0, 0, bs);

% element center and size
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);

nv = length(u0);
U0 = zeros(bs.nb, nv * msh.nLElems);
switch layout   
    case 1
        for i = 1 : nv
            U0(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems) = ME \ (bs.phitw{1} * u0{i}(ct + 0.5 * h .* quad.points));
        end
    case 2
        for i = 1 : nv
            U0(:, i : nv : end) = ME \ (bs.phitw{1} * u0{i}(ct + 0.5 * h .* quad.points));
        end
end

end

