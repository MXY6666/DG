% trial_der : derivative order of trial basis functions
% test_der  : derivative order of test basis functions
% bs        : basis function set data (a struct) or just the degree of 
%             polynomial (an integer)
% mtol      : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function EM = computeElementMatrix_refCube(trial_der, test_der, bs, mtol)

if (nargin < 1) || isempty(trial_der)
    trial_der = [0, 0, 0];
end

if (nargin < 2) || isempty(test_der)
    test_der = [0, 0, 0];
end

if (length(trial_der) ~= 3) || (length(test_der) ~= 3)
    error('Wrong size of derivative order')
end

if all(trial_der(:) == [0; 0; 0])
    trial_ind = 1;
elseif all(trial_der(:) == [1; 0; 0])
    trial_ind = 2;
elseif all(trial_der(:) == [0; 1; 0])
    trial_ind = 3;
elseif all(trial_der(:) == [0; 0; 1])
    trial_ind = 4;    
else
    error('Not implemented derivative order of trial basis funtion')
end

if all(test_der(:) == [0; 0; 0])
    test_ind = 1;
elseif all(test_der(:) == [1; 0; 0])
    test_ind = 2;
elseif all(test_der(:) == [0; 1; 0])
    test_ind = 3;
elseif all(test_der(:) == [0; 0; 1])
    test_ind = 4;     
else
    error('Not implemented derivative order of trial basis funtion')
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
    quad       = GaussQuadratureRule_cube([k + 1, k + 1, k + 1], 301);
    phi_trial  = basisFunctionSet_cube(quad.points(:, 1), quad.points(:, 2), quad.points(:, 3), 300 + k, trial_der);
    phi_test   = basisFunctionSet_cube(quad.points(:, 1), quad.points(:, 2), quad.points(:, 3), 300 + k, test_der);  
    phitw_test = (quad.weights .* phi_test)';  
    EM         = phitw_test * phi_trial;
elseif isstruct(bs)
    if ~strcmpi(bs.refGeom, 'cube')
        error('Wrong reference geometry for basis functions to evaluate on')
    end
    if (bs.elemPointsType == 301 && 2 * bs.neps(1) - 1 >= 2 * bs.deg && 2 * bs.neps(2) - 1 >= 2 * bs.deg && 2 * bs.neps(3) - 1 >= 2 * bs.deg) 
        EM = bs.phitw{test_ind} * bs.phi{trial_ind};
    else      
        quad       = GaussQuadratureRule_cube([bs.deg + 1, bs.deg + 1, bs.deg + 1], 301);
        phi_trial  = basisFunctionSet_cube(quad.points(:, 1), quad.points(:, 2), quad.points(:, 3), bs.type, trial_der);
        phi_test   = basisFunctionSet_cube(quad.points(:, 1), quad.points(:, 2), quad.points(:, 3), bs.type, test_der);  
        phitw_test = (quad.weights .* phi_test)';  
        EM         = phitw_test * phi_trial;
    end    
else 
    error('Wrong argument bs')     
end
EM = mychop(EM, mtol);

end




