% trial_der : derivative order of trial basis functions
% test_der  : derivative order of test basis functions
% bs        : basis function set data (a struct) or just the degree of 
%             polynomial (an integer)
% mtol      : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function EM = computeElementMatrix_refLine(trial_der, test_der, bs, mtol)

if (nargin < 1) || isempty(trial_der)
    trial_der = 0;
end

if (nargin < 2) || isempty(test_der)
    test_der = 0;
end

if (trial_der ~= 0) && (trial_der ~= 1)
    error('Not implemented derivative order of trial basis funtion')
end

if (test_der ~= 0) && (test_der ~= 1)
    error('Not implemented derivative order of test basis funtion')
end

% Default we choose basis functions of degree 1 
if (nargin < 3) || isempty(bs)
    bs = 1;
end

if (nargin < 4) || isempty(mtol)
    mtol = 1.0e-12;
end

if isfloat(bs) && (bs >= 0)
    k          = floor(bs);
    quad       = GaussQuadratureRule_line(k + 1, 101);
    phi_trial  = basisFunctionSet_line(quad.points, 100 + k, trial_der);
    phi_test   = basisFunctionSet_line(quad.points, 100 + k, test_der);
    phitw_test = (quad.weights .* phi_test)';
    EM         = phitw_test * phi_trial;
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'Line')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if (bs.elemPointsType == 101 && 2 * bs.nep - 1 >= 2 * bs.deg) || ...
       (bs.elemPointsType == 102 && 2 * bs.nep - 3 >= 2 * bs.deg)
        EM = bs.phitw{test_der + 1} * bs.phi{trial_der + 1};
    else
        quad       = GaussQuadratureRule_line(bs.deg + 1, 101);
        phi_trial  = basisFunctionSet_line(quad.points, bs.type, trial_der);
        phi_test   = basisFunctionSet_line(quad.points, bs.type, test_der);
        phitw_test = (quad.weights .* phi_test)';
        EM         = phitw_test * phi_trial;
    end        
else 
    error('Wrong argument bs')    
end
EM = mychop(EM, mtol);

end


