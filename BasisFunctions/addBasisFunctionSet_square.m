% Evaluate basis functions at Gauss points of four subelements of reference 
% square and eight subfaces of reference element boundaries
%               F11    F12
%            o------o------o
%            |      |      |
%         F6 |  3   |   4  | F8
%            |      |      |
%            o------o------o
%            |      |      |
%         F5 |  1   |   2  | F7
%            |      |      |
%            o------o------o
%               F9     F10
% Copyright (C) Xiangyi Meng.

function bs = addBasisFunctionSet_square(bs, quad1, quad2)

bs.phi_sube = {basisFunctionSet_square(-0.5 + 0.5 * quad2.points(:, 1), -0.5 + 0.5 * quad2.points(:, 2), bs.type, [0, 0]), ...
               basisFunctionSet_square( 0.5 + 0.5 * quad2.points(:, 1), -0.5 + 0.5 * quad2.points(:, 2), bs.type, [0, 0]), ...
               basisFunctionSet_square(-0.5 + 0.5 * quad2.points(:, 1),  0.5 + 0.5 * quad2.points(:, 2), bs.type, [0, 0]), ...
               basisFunctionSet_square( 0.5 + 0.5 * quad2.points(:, 1),  0.5 + 0.5 * quad2.points(:, 2), bs.type, [0, 0])};
           
bs.phitw_sube = cell(1, 4);
for j = 1 : 4
    bs.phitw_sube{j} = (quad2.weights .* bs.phi_sube{j})';
end         

if (bs.faceOrderType == 1)             
    bs.phi_face(1, 5 : 12) = {basisFunctionSet_square(-1, -0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...
                              basisFunctionSet_square(-1,  0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 1, -0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 1,  0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...
                              basisFunctionSet_square(-0.5 + 0.5 * quad1.points, -1, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 0.5 + 0.5 * quad1.points, -1, bs.type, [0, 0]), ...
                              basisFunctionSet_square(-0.5 + 0.5 * quad1.points,  1, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 0.5 + 0.5 * quad1.points,  1, bs.type, [0, 0])};
elseif (bs.faceOrderType == 2)
    bs.phi_face(1, 5 : 12) = {basisFunctionSet_square(-1,  0.5 + 0.5 * flip(quad1.points), bs.type, [0, 0]), ...
                              basisFunctionSet_square(-1, -0.5 + 0.5 * flip(quad1.points), bs.type, [0, 0]), ...
                              basisFunctionSet_square( 1, -0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 1,  0.5 + 0.5 * quad1.points, bs.type, [0, 0]), ...                              
                              basisFunctionSet_square(-0.5 + 0.5 * quad1.points, -1, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 0.5 + 0.5 * quad1.points, -1, bs.type, [0, 0]), ...
                              basisFunctionSet_square( 0.5 + 0.5 * flip(quad1.points),  1, bs.type, [0, 0]), ...
                              basisFunctionSet_square(-0.5 + 0.5 * flip(quad1.points),  1, bs.type, [0, 0])};
else
    error('Wrong face order type')
end
for j = 5 : 12
    bs.phi_face{2, j} = flip(bs.phi_face{1, j});
end

bs.phitw_face(:, 5 : 12) = cell(2, 8);
for j = 5 : 12
    bs.phitw_face{1, j} = (quad1.weights .* bs.phi_face{1, j})';
    bs.phitw_face{2, j} = (quad1.weights .* bs.phi_face{2, j})';
end


end


