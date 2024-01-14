% Copyright (C) Xiangyi Meng.

function bs = setBasisFunctionSet_cube(quad2, quad3, basisType)

% reference geometry to evaluate
bs.refGeom = 'cube';

% type of basis functions
bs.type = basisType;

% the maximum polynomial degree w.r.t. one variable
bs.deg = mod(basisType, 10);

% nodes tag for evaluation
bs.elemPointsType = quad3.type;
bs.facePointsType = quad2.type;

% Evalute basis functions and their derivatives at element Gauss nodes of 
% reference element [-1, 1]^3
bs.phi = {basisFunctionSet_cube(quad3.points(:, 1), quad3.points(:, 2), quad3.points(:, 3), basisType, [0, 0, 0]), ...
          basisFunctionSet_cube(quad3.points(:, 1), quad3.points(:, 2), quad3.points(:, 3), basisType, [1, 0, 0]), ...
          basisFunctionSet_cube(quad3.points(:, 1), quad3.points(:, 2), quad3.points(:, 3), basisType, [0, 1, 0]), ...
          basisFunctionSet_cube(quad3.points(:, 1), quad3.points(:, 2), quad3.points(:, 3), basisType, [0, 0, 1])};
      
% Evalute basis functions at face Gauss nodes of the refence element 
% [-1, 1]^3
bs.phi_face = cell(8, 6);
for i = 1 : 8
    switch i 
        case 1
            coor1 = quad2.points(:, 1);
            coor2 = quad2.points(:, 2);
        case 2
            coor1 = quad2.points(:, 2);
            coor2 = -quad2.points(:, 1);
        case 3
            coor1 = -quad2.points(:, 1);
            coor2 = -quad2.points(:, 2);            
        case 4
            coor1 = -quad2.points(:, 2);
            coor2 = quad2.points(:, 1);            
        case 5
            coor1 = -quad2.points(:, 1);
            coor2 = quad2.points(:, 2);            
        case 6
            coor1 = quad2.points(:, 1);
            coor2 = -quad2.points(:, 2);            
        case 7
            coor1 = quad2.points(:, 2);
            coor2 = quad2.points(:, 1);            
        case 8
            coor1 = -quad2.points(:, 2);
            coor2 = -quad2.points(:, 1);            
    end
    bs.phi_face(i, :) = {basisFunctionSet_cube(-1, coor1, coor2, basisType, [0, 0, 0]), ...
                         basisFunctionSet_cube( 1, coor1, coor2, basisType, [0, 0, 0]), ...
                         basisFunctionSet_cube(coor2, -1, coor1, basisType, [0, 0, 0]), ...
                         basisFunctionSet_cube(coor2,  1, coor1, basisType, [0, 0, 0]), ...
                         basisFunctionSet_cube(coor1, coor2, -1, basisType, [0, 0, 0]), ...
                         basisFunctionSet_cube(coor1, coor2,  1, basisType, [0, 0, 0])};       
end
                 
% Weights multiply basis functions and their derivatives componentwisely to 
% prepare for element integration
bs.phitw = {(quad3.weights .* bs.phi{1})', ...
            (quad3.weights .* bs.phi{2})', ...
            (quad3.weights .* bs.phi{3})', ...
            (quad3.weights .* bs.phi{4})'};
         
% Weights multiply basis functions componentwisely to prepare for face 
% integration
bs.phitw_face = cell(8, 6);
for i = 1 : 8
    for j = 1 : 6
       bs.phitw_face{i, j} = (quad2.weights .* bs.phi_face{i, j})';
    end
end   

% number of basis functions and its square
bs.nb = size(bs.phi{1}, 2);
bs.nb2 = bs.nb * bs.nb;

% number of Gauss points at element and face
bs.nep = quad3.np;
bs.nfp = quad2.np;
bs.neps = quad3.nps;
bs.nfps = quad2.nps;
        
end



