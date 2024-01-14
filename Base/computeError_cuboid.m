% msh    : mesh of cuboid element in 3D
% u      : cell array of exact solutions
% U      : coefficients of numerical solutions
% t      : time
% bs     : basis function set data
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = computeError_cuboid(msh, u, U, t, bs, layout)

if (nargin < 5)
    error('Not enough arguments')
end

if (msh.type ~= 301) && (msh.type ~= 302)
    error('Wrong mesh type')
end

if (size(U, 1) ~= bs.nb)
    error('The number of rows of U does not match')
end

nv = length(u);  
if (size(U, 2) ~= nv * msh.nLElems)
    error('The number of columns of U does not match')
end

if (nargin < 6) || isempty(layout)
    layout = 1;
end
if (layout ~= 1) && (layout ~= 2)
    error('Wrong argument layout')
end

% Generate quadrature rule based on bs
quad3 = GaussQuadratureRule_cube(bs.neps, 301); 
if ~(bs.elemPointsType == 301 && bs.neps(1) >= bs.deg + 1 && bs.neps(2) >= bs.deg + 1 && bs.neps(3) >= bs.deg + 1)
    quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
    quad3 = GaussQuadratureRule_cube([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
    bs = setBasisFunctionSet_cube(quad2, quad3, bs.type);
end 

% center, size and Jacobian in each leaf element
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);
J  = msh.elemJac(:, msh.LElems);

L2err = zeros(nv, 1);
Lierr = zeros(nv, 1);
for i = 1 : nv
    % Evaluate exact solutions on Gauss quadrature nodes
    exactSolutions = u{i}(ct(1, :) + 0.5 * h(1, :) .* quad3.points(:, 1), ct(2, :) + 0.5 * h(2, :) .* quad3.points(:, 2), ct(3, :) + 0.5 * h(3, :) .* quad3.points(:, 3), t);
    
    switch layout
        case 1
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems);            
        case 2
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, i : nv : end);            
    end
    
    % Compute L2 error
    L2err(i) = sqrt(sum(J .* (quad3.weights' * (exactSolutions - numericalSolutions).^2)));  

    % Compute L_infinity error
    Lierr(i) = max(max(abs(exactSolutions - numericalSolutions)));
end


end



