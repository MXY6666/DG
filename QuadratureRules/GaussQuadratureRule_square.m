% Copyright (C) Xiangyi Meng.

function quad = GaussQuadratureRule_square(np, type)

if (numel(np) ~= 2)
    error('The first argument should be a 2D vector')
end

quad.type = type;
quad.np = prod(np); quad.nps = np;

% type 201 : Gauss Legendre in both x and y direction
% type 202 : Gauss Lobatto in both x and y direction
% type 203 : Gauss Legendre and Gauss Lobatto in x and y direction, respectively
% type 204 : Gauss Lobatto and Gauss Legendre in x and y direction, respectively
switch type
    case 201
        quad1 = GaussQuadratureRule_line(np(1), 101);
        quad2 = GaussQuadratureRule_line(np(2), 101);
    case 202
        quad1 = GaussQuadratureRule_line(np(1), 102);
        quad2 = GaussQuadratureRule_line(np(2), 102);
    case 203
        quad1 = GaussQuadratureRule_line(np(1), 101);
        quad2 = GaussQuadratureRule_line(np(2), 102);        
    case 204
        quad1 = GaussQuadratureRule_line(np(1), 102);
        quad2 = GaussQuadratureRule_line(np(2), 101);          
    otherwise
        error('Not implemented type of Gauss quadrature rule')
end
quad.points = [repmat(quad1.points, [np(2), 1]), repelem(quad2.points, np(1))];
quad.weights = reshape(quad1.weights * quad2.weights', [quad.np, 1]);


end