% msh    : mesh of rectangular element in 2D
% u      : cell array of exact solutions
% U      : coefficients of numerical solutions
% t      : time
% bs     : basis function set data
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = computeError_rect(msh, u, U, t, bs, layout)

if (nargin < 5)
    error('Not enough arguments')
end

if (msh.type ~= 201) && (msh.type ~= 202)
    error('Wrong mesh type')
end

if (size(U, 1) ~= bs.nb)
    error('The number of rows of U does not match')
end

nv = length(u);  
if (size(U, 2) ~= nv * msh.nLElems)
    error('The number of columns of U does not match')
end

nb = size(U, 1); % number of basis functions
deg1 = (sqrt(8 * nb + 1) - 3) / 2;
deg2 = sqrt(nb) - 1;
if (mod(deg1, 1) == 0)
    if (bs.type ~= 200 + deg1)
        error('The arguments U and bs are not consistent')
    end    
elseif (mod(deg2, 1) == 0)   
    if (bs.type ~= 210 + deg2 && bs.type ~= 220 + deg2)
        error('The arguments U and bs are not consistent')
    end       
else
    error('Wrong number of rows of argument U')
end

if (nargin < 6) || isempty(layout)
    layout = 1;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

% Generate quadrature rule based on bs
quad2 = GaussQuadratureRule_square(bs.neps, bs.elemPointsType);
if (bs.elemPointsType == 201) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 1)
    quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
    quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
    bs = setBasisFunctionSet_square(quad1, quad2, bs.type);
end
if (bs.elemPointsType == 202) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 2)
    quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
    quad2 = GaussQuadratureRule_square([bs.deg + 2, bs.deg + 2], 202);
    bs = setBasisFunctionSet_square(quad1, quad2, bs.type);
end       
if (bs.elemPointsType == 203) && (bs.neps(1) < bs.deg + 1 || bs.neps(2) < bs.deg + 2)
    quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
    quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 2], 203);
    bs = setBasisFunctionSet_square(quad1, quad2, bs.type);
end     
if (bs.elemPointsType == 204) && (bs.neps(1) < bs.deg + 2 || bs.neps(2) < bs.deg + 1)
    quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
    quad2 = GaussQuadratureRule_square([bs.deg + 2, bs.deg + 1], 204);
    bs = setBasisFunctionSet_square(quad1, quad2, bs.type);
end

% center, size and Jacobian in each leaf element
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);
J  = msh.elemJac(:, msh.LElems);

L2err = zeros(nv, 1);
Lierr = zeros(nv, 1);
for i = 1 : nv
    % Evaluate exact solutions on Gauss quadrature nodes
    exactSolutions = u{i}(ct(1, :) + 0.5 * h(1, :) .* quad2.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad2.points(:, 2), t);
    
    switch layout
        case 1
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems);            
        case 2
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, i : nv : end);            
    end
    
    % Compute L2 error
    L2err(i) = sqrt(sum(J .* (quad2.weights' * (exactSolutions - numericalSolutions).^2)));  

    % Compute L_infinity error
    Lierr(i) = max(max(abs(exactSolutions - numericalSolutions)));
end

end

