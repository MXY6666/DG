% Copyright (C) Xiangyi Meng.

function quad = GaussQuadratureRule_cube(np, type)

if (numel(np) ~= 3)
    error('The first argument should be a 3D vector')
end

quad.type = type;
quad.np = prod(np); quad.nps = np;

% type 301 : Gauss Legendre in all the direction
switch type
    case 301
        quad1 = GaussQuadratureRule_line(np(1), 101);
        quad2 = GaussQuadratureRule_line(np(2), 101);
        quad3 = GaussQuadratureRule_line(np(3), 101);
    otherwise
        error('Not implemented type of Gauss quadrature rule')         
end
quad.points = [repmat(quad1.points, [np(2) * np(3), 1]), repmat(repelem(quad2.points, np(1)), [np(3), 1]), repelem(quad3.points, np(1) * np(2))];
quad.weights = reshape(reshape(quad1.weights * quad2.weights', [np(1) * np(2), 1]) * quad3.weights', [quad.np, 1]);

end