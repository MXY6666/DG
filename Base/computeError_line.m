% msh    : mesh of line element in 1D
% u      : cell array of exact solutions
% U      : coefficients of numerical solutions
% t      : time
% bs     : basis function set data
% layout : 1 or 2, layout of U0
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = computeError_line(msh, u, U, t, bs, layout)

if (nargin < 5)
    error('Not enough arguments')
end

if (msh.type ~= 101) && (msh.type ~= 102)
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
quad = GaussQuadratureRule_line(bs.nep, bs.elemPointsType);
if (bs.elemPointsType == 101) && (bs.nep < bs.deg + 1)                                                     
    quad = GaussQuadratureRule_line(bs.deg + 1, 101);
    bs   = setBasisFunctionSet_line(quad, bs.type);
end
if (bs.elemPointsType == 102) && (bs.nep < bs.deg + 2)
    quad = GaussQuadratureRule_line(bs.deg + 2, 102);
    bs   = setBasisFunctionSet_line(quad, bs.type);                                                                            
end

% center, size and Jacobian in each leaf element
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);
J  = msh.elemJac(:, msh.LElems);

nv = length(u);  
if (size(U, 2) ~= nv * msh.nLElems)
    error('The number of columns of U does not match')
end
L2err = zeros(nv, 1);
Lierr = zeros(nv, 1);
for i = 1 : nv
    % Eval exact solutions on Gauss quadrature nodes
    exactSolutions = u{i}(ct + 0.5 * h .* quad.points, t);
    
    switch layout
        case 1
            % Eval numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, (i - 1) * msh.nLElems + 1 : i * msh.nLElems);            
        case 2
            % Eval numerical solutions on Gauss quadrature nodes
            numericalSolutions = bs.phi{1} * U(:, i : nv : end);            
    end
    
    % Compute L2 error
    L2err(i) = sqrt(sum(J .* (quad.weights' * (exactSolutions - numericalSolutions).^2))); 

    % Compute L_infinity error
    Lierr(i) = max(max(abs(exactSolutions - numericalSolutions)));   
end
 
end

