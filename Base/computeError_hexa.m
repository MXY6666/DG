% msh    : mesh of linear hexahedral element in 2D
% u      : cell array of exact solutions
% U      : coefficients of numerical solutions
% t      : time
% bs     : basis function set data
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = computeError_hexa(msh, u, U, t, bs, layout)

if (nargin < 5)
    error('Not enough arguments')
end

if (msh.type ~= 311) && (msh.type ~= 312)
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
quad3 = GaussQuadratureRule_hexa(bs.neps, 301); 
if ~(bs.elemPointsType == 301 && bs.neps(1) >= bs.deg + 1 && bs.neps(2) >= bs.deg + 1 && bs.neps(3) >= bs.deg + 1)
    quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
    quad3 = GaussQuadratureRule_hexa([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
    bs = setBasisFunctionSet_cube(quad2, quad3, bs.type);
end 

if ~all(msh.nGPs == bs.neps)
    error('The number of Gauss ponits in msh and bs does not coinside')
end

% Gauss points and Jacobian in each leaf element
gpx = msh.elemGPx(:, msh.LElems);
gpy = msh.elemGPy(:, msh.LElems);
gpz = msh.elemGPz(:, msh.LElems);
J   = msh.elemJac(:, msh.LElems);

L2err = zeros(nv, 1);
Lierr = zeros(nv, 1);
for i = 1 : nv
    % Evaluate exact solutions on Gauss quadrature nodes
    exactSolutions = u{i}(gpx, gpy, gpz, t);
    
    switch layout
        case 1
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems);            
        case 2
            % Evaluate numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, i : nv : end);            
    end
    
    % Compute L2 error
    L2err(i) = sqrt(sum(quad3.weights' * ((exactSolutions - numericalSolutions).^2 .* J)));  

    % Compute L_infinity error
    Lierr(i) = max(max(abs(exactSolutions - numericalSolutions)));
end


end




