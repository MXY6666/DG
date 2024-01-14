% The isoparametric mapping form reference square to a quadrilateral with
% vertices (x_i, y_i), i = 1, 2, 3, 4 is
% (x, y) = \sum{i = 1}^{4}(x_i, y_i)N_i(\xi, \eta), 
% where
% N_1 = (1 - \xi)(1 - \eta) / 4, 
% N_2 = (1 + \xi)(1 - \eta) / 4,
% N_3 = (1 + \xi)(1 + \eta) / 4, 
% N_4 = (1 - \xi)(1 + \eta) / 4.
% Note, the points stored in quad are actually not the Gauss points in the
% physical quadrilateral, but we can get it through 
% quad.N * (x1, y1; x2, y2; x3, y3; x4, y4) as long as the coordinates of
% the quadrilateral vertices are provided
% Copyright (C) Xiangyi Meng.

function quad = GaussQuadratureRule_quad(np, type)

quad = GaussQuadratureRule_square(np, type);

quad.N = 0.25 * [(1 - quad.points(:, 1)) .* (1 - quad.points(:, 2)), ...
                 (1 + quad.points(:, 1)) .* (1 - quad.points(:, 2)), ...
                 (1 + quad.points(:, 1)) .* (1 + quad.points(:, 2)), ...
                 (1 - quad.points(:, 1)) .* (1 + quad.points(:, 2))];
             
quad.N_xi = 0.25 * [quad.points(:, 2) - 1, 1 - quad.points(:, 2), ...
                    1 + quad.points(:, 2), -(1 + quad.points(:, 2))];
                
quad.N_eta = 0.25 * [quad.points(:, 1) - 1, -(1 + quad.points(:, 1)), ...
                     1 + quad.points(:, 1), 1 - quad.points(:, 1)];

end