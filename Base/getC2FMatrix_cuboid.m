% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function c2f = getC2FMatrix_cuboid(bs, mtol)

% Default we choose basis functions of degree 1 
if (nargin < 1) || isempty(bs)
    bs = 1;
end

if (nargin < 2) || isempty(mtol)
    mtol = 1.0e-12;
end

if isfloat(bs) && (bs >= 1)
    k = floor(bs);
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    quad3 = GaussQuadratureRule_cube([k + 1, k + 1, k + 1], 301);
    bs    = setBasisFunctionSet_cube(quad2, quad3, 300 + k);
    bs    = addBasisFunctionSet_cube(bs, quad2, quad3);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'cube')
        error('Wrong reference geometry for basis functions to evaluate on')
    end   
    if ~(bs.elemPointsType == 301 && bs.neps(1) >= bs.deg + 1 && bs.neps(2) >= bs.deg + 1 && bs.neps(3) >= bs.deg + 1)
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        quad3 = GaussQuadratureRule_cube([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
        bs    = setBasisFunctionSet_cube(quad2, quad3, bs.type);
        bs    = addBasisFunctionSet_cube(bs, quad2, quad3);
    end  
else
    error('Wrong argument bs')
end

m  = bs.phitw{1} * bs.phi{1};
dm = decomposition(m);

c2f = cell(1, 8);
for i = 1 : 8
    c2f{i} = mychop(dm \ (bs.phitw{1} * bs.phi_sube{i}), mtol);
end

end

