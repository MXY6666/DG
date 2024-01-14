% Evaluate basis functions at Gauss points of eight subelements of 
% reference cube and 24 subfaces of reference cube boundaries
% Copyright (C) Xiangyi Meng.

function bs = addBasisFunctionSet_cube(bs, quad2, quad3)

bs.phi_sube = {basisFunctionSet_cube(-0.5 + 0.5 * quad3.points(:, 1), -0.5 + 0.5 * quad3.points(:, 2), -0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube( 0.5 + 0.5 * quad3.points(:, 1), -0.5 + 0.5 * quad3.points(:, 2), -0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube(-0.5 + 0.5 * quad3.points(:, 1),  0.5 + 0.5 * quad3.points(:, 2), -0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube( 0.5 + 0.5 * quad3.points(:, 1),  0.5 + 0.5 * quad3.points(:, 2), -0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...               
               basisFunctionSet_cube(-0.5 + 0.5 * quad3.points(:, 1), -0.5 + 0.5 * quad3.points(:, 2),  0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube( 0.5 + 0.5 * quad3.points(:, 1), -0.5 + 0.5 * quad3.points(:, 2),  0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube(-0.5 + 0.5 * quad3.points(:, 1),  0.5 + 0.5 * quad3.points(:, 2),  0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0]), ...
               basisFunctionSet_cube( 0.5 + 0.5 * quad3.points(:, 1),  0.5 + 0.5 * quad3.points(:, 2),  0.5 + 0.5 * quad3.points(:, 3), bs.type, [0, 0, 0])};
bs.phitw_sube = cell(1, 8);
for i = 1 : 8
    bs.phitw_sube{i} = (quad3.weights .* bs.phi_sube{i})';
end

bs.phi_face(:, 7 : 30) = cell(8, 24);
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
    bs.phi_face{i, 7}  = basisFunctionSet_cube(-1, -0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 8}  = basisFunctionSet_cube(-1,  0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 9}  = basisFunctionSet_cube(-1, -0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 10} = basisFunctionSet_cube(-1,  0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
        
    bs.phi_face{i, 11} = basisFunctionSet_cube( 1, -0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 12} = basisFunctionSet_cube( 1,  0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 13} = basisFunctionSet_cube( 1, -0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
    bs.phi_face{i, 14} = basisFunctionSet_cube( 1,  0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, bs.type, [0, 0, 0]);
            
    bs.phi_face{i, 15} = basisFunctionSet_cube(-0.5 + 0.5 * coor2, -1, -0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 16} = basisFunctionSet_cube(-0.5 + 0.5 * coor2, -1,  0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 17} = basisFunctionSet_cube( 0.5 + 0.5 * coor2, -1, -0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 18} = basisFunctionSet_cube( 0.5 + 0.5 * coor2, -1,  0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    
    bs.phi_face{i, 19} = basisFunctionSet_cube(-0.5 + 0.5 * coor2,  1, -0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 20} = basisFunctionSet_cube(-0.5 + 0.5 * coor2,  1,  0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 21} = basisFunctionSet_cube( 0.5 + 0.5 * coor2,  1, -0.5 + 0.5 * coor1, bs.type, [0, 0, 0]); 
    bs.phi_face{i, 22} = basisFunctionSet_cube( 0.5 + 0.5 * coor2,  1,  0.5 + 0.5 * coor1, bs.type, [0, 0, 0]);
      
    bs.phi_face{i, 23} = basisFunctionSet_cube(-0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, -1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 24} = basisFunctionSet_cube( 0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2, -1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 25} = basisFunctionSet_cube(-0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, -1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 26} = basisFunctionSet_cube( 0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2, -1, bs.type, [0, 0, 0]);
       
    bs.phi_face{i, 27} = basisFunctionSet_cube(-0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2,  1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 28} = basisFunctionSet_cube( 0.5 + 0.5 * coor1, -0.5 + 0.5 * coor2,  1, bs.type, [0, 0, 0]);
    bs.phi_face{i, 29} = basisFunctionSet_cube(-0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2,  1, bs.type, [0, 0, 0]); 
    bs.phi_face{i, 30} = basisFunctionSet_cube( 0.5 + 0.5 * coor1,  0.5 + 0.5 * coor2,  1, bs.type, [0, 0, 0]);      
end

bs.phitw_face(:, 7 : 30) = cell(8, 24);
for i = 1 : 8
    for j = 7 : 30
       bs.phitw_face{i, j} = (quad2.weights .* bs.phi_face{i, j})';
    end
end 
             
end