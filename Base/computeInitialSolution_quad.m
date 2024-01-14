% msh    : mesh of linear quadrilateral element in 2D
% u0     : cell array of exact initial solutions
% bs     : basis function set data (a struct) or just the degree of 
%          polynomial (an integer)
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function U0 = computeInitialSolution_quad(msh, u0, bs, layout)

if (nargin < 2)
    error('Not enough arguments')
end

if (msh.type ~= 211) && (msh.type ~= 212)
    error('Wrong mesh type')
end

if (nargin < 3) || isempty(bs)
    bs = 1;
end

if isfloat(bs) && (bs >= 1)
    k     = floor(bs);
    quad1 = GaussQuadratureRule_line(k + 1, 101); 
    quad2 = GaussQuadratureRule_quad(msh.nGPs, 201);
    bs    = setBasisFunctionSet_square(quad1, quad2, 200 + k, 2);
end
if isstruct(bs)
    if ~strcmpi(bs.refGeom, 'square')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if (bs.elemPointsType == 201) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 1)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_quad([bs.deg + 1, bs.deg + 1], 201);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type, 2);
    end
    if (bs.elemPointsType == 202) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 2)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_quad([bs.deg + 2, bs.deg + 2], 202);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type, 2);
    end       
    if (bs.elemPointsType == 203) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 2)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_quad([bs.deg + 1, bs.deg + 2], 203);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type, 2);
    end     
    if (bs.elemPointsType == 204) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 1)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_quad([bs.deg + 2, bs.deg + 1], 204);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type, 2);
    end 
else
    error('Wrong argument bs')
end

if ~all(msh.nGPs == bs.neps) || (size(msh.massMatInv{1}, 1) ~= bs.nb)
    error('msh and bs does not coinside')
end

if (nargin < 4) || isempty(layout)
    layout = 1;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

% Gauss points and Jacobian in each leaf element
gpx = msh.elemGPx(:, msh.LElems);
gpy = msh.elemGPy(:, msh.LElems);
J   = msh.elemJac(:, msh.LElems);

nv = length(u0);
U0 = zeros(bs.nb, nv * msh.nLElems);
switch layout
    case 1
        for i = 1 : nv
            U0(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems) = cell2mat(cellfun(@mmult, msh.massMatInv(:, msh.LElems), mat2cell(bs.phitw{1} * (u0{i}(gpx, gpy) .* J), bs.nb, ones(1, msh.nLElems)), 'UniformOutput', false));            
        end    
    case 2
        for i = 1 : nv
            U0(:, i : nv : end) = cell2mat(cellfun(@mmult, msh.massMatInv(:, msh.LElems), mat2cell(bs.phitw{1} * (u0{i}(gpx, gpy) .* J), bs.nb, ones(1, msh.nLElems)), 'UniformOutput', false));            
        end 
end


end


