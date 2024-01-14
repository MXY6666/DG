% mtol : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function f2c = getF2CMatrix_line(bs, mtol)

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
    bs    = setBasisFunctionSet_line(quad1, 100 + k);
    bs    = addBasisFunctionSet_line(bs, quad1);
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'line')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if ~(bs.elemPointsType == 101 && 2 * bs.nep - 1 >= 2 * bs.deg)
        quad1 = GaussQuadratureRule_line(bs.deg + 1, 101);
        bs    = setBasisFunctionSet_line(quad1, bs.type);
        bs    = addBasisFunctionSet_line(bs, quad1);
    end
else
     error('Wrong argument bs')
end

m  = bs.phitw{1} * bs.phi{1};
dm = decomposition(m);

f2c = cell(1, 2);
for i = 1 : 2
    f2c{i} = 0.5 * mychop(dm \ (bs.phitw_sube{i} * bs.phi{1}), mtol);
end

end

