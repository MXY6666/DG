% facePointsOrder = 1 : the four faces are ordered by left, right, bottom and 
% top, and the points on the four face are placed in an increasing order 
% w.r.t. x or y
% facePointsOrder = 2 : the four faces are ordered by left, right, bottom and
% top, and the points on the four face are placed in a counterclockwise 
% order
% Copyright (C) Xiangyi Meng.

function bs = setBasisFunctionSet_square(quad1, quad2, basisType, facePointsOrder)

if (nargin < 3)
    error('Not enough arguments')
end

if (nargin < 4) || isempty(facePointsOrder)
    facePointsOrder = 1;
end
if (facePointsOrder ~= 1) && (facePointsOrder ~= 2)
    error('Wrong face points order')
end 

% reference geometry to evaluate
bs.refGeom = 'square';

% type of basis functions
bs.type = basisType;

% the maximum polynomial degree w.r.t. one variable
bs.deg = mod(basisType, 10);

% nodes tag for evaluation
bs.elemPointsType = quad2.type;
bs.facePointsType = quad1.type;

% type of face points order
bs.faceOrderType = facePointsOrder;

% Evalute basis functions and their derivatives at element Gauss nodes of 
% reference element [-1, 1]^2
bs.phi = {basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [0, 0]), ...
          basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [1, 0]), ...
          basisFunctionSet_square(quad2.points(:, 1), quad2.points(:, 2), basisType, [0, 1])};

% Evalute basis functions at face Gauss nodes of the refence element 
% [-1, 1]^2
bs.phi_face = cell(2, 4);
if (facePointsOrder == 1)
    bs.phi_face(1, :) = {basisFunctionSet_square(-1, quad1.points, basisType, [0, 0]), ...
                         basisFunctionSet_square( 1, quad1.points, basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points, -1, basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points,  1, basisType, [0, 0])};    
else
    bs.phi_face(1, :) = {basisFunctionSet_square(-1, flip(quad1.points), basisType, [0, 0]), ...
                         basisFunctionSet_square( 1, quad1.points,       basisType, [0, 0]), ...
                         basisFunctionSet_square(quad1.points,      -1,  basisType, [0, 0]), ...
                         basisFunctionSet_square(flip(quad1.points), 1,  basisType, [0, 0])};    
end           
bs.phi_face(2, :) = {flip(bs.phi_face{1, 1}), flip(bs.phi_face{1, 2}), flip(bs.phi_face{1, 3}), flip(bs.phi_face{1, 4})}; 

% Weights multiply basis functions and their derivatives componentwisely to 
% prepare for element integration
bs.phitw = {(quad2.weights .* bs.phi{1})', ...
            (quad2.weights .* bs.phi{2})', ...
            (quad2.weights .* bs.phi{3})'};

% Weights multiply basis functions componentwisely to prepare for face 
% integration
bs.phitw_face = cell(2, 4);
for i = 1 : 2
    for j = 1 : 4
       bs.phitw_face{i, j} = (quad1.weights .* bs.phi_face{i, j})';
    end
end

% number of basis functions and its square
bs.nb = size(bs.phi{1}, 2);
bs.nb2 = bs.nb * bs.nb;

% number of Gauss points at element and face
bs.nep = quad2.np;
bs.nfp = quad1.np;
bs.neps = quad2.nps;

end

