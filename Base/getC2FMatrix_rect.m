% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function c2f = getC2FMatrix_rect(bs, mtol)

% Default we choose basis functions of degree 1 
if (nargin < 1) || isempty(bs)
    bs = 1;
end

if (nargin < 2) || isempty(mtol)
    mtol = 1.0e-12;
end

if isfloat(bs) && (bs >= 1)
    k     = floor(bs);
    quad1 = GaussQuadratureRule_line(k + 1, 101);
    quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);
    bs    = setBasisFunctionSet_square(quad1, quad2, 200 + k);
    bs    = addBasisFunctionSet_square(bs, quad1, quad2);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'square')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if ~(bs.elemPointsType == 201 && 2 * bs.neps(1) - 1 >= 2 * bs.deg && 2 * bs.neps(2) - 1 >= 2 * bs.deg)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        quad2 = GaussQuadratureRule_square([bs.deg + 1, bs.deg + 1], 201);
        bs    = setBasisFunctionSet_square(quad1, quad2, bs.type);
        bs    = addBasisFunctionSet_square(bs, quad1, quad2);
    end
else
     error('Wrong argument bs')
end

m  = bs.phitw{1} * bs.phi{1};
dm = decomposition(m);

c2f = cell(1, 4);
for i = 1 : 4
    c2f{i} = mychop(dm \ (bs.phitw{1} * bs.phi_sube{i}), mtol);
end


end


