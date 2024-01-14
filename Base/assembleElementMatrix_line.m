% msh       : mesh of line element in 1D
% trial_der : derivative order of trial basis functions
% test_der  : derivative order of test basis functions
% bs        : basis function set data (a struct) or just the degree of 
%             polynomial (an integer)
% inverse   : inverse the element matrix or not 
% mtol      : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function A = assembleElementMatrix_line(msh, trial_der, test_der, bs, inverse, EM, mtol)

if (nargin < 4)
    error('Not enough arguments')
end    

if (msh.type ~= 101) && (msh.type ~= 102)
    error('Wrong mesh type')
end

if (nargin < 5) || isempty(inverse)
    inverse = false;
end

if (nargin < 7) || isempty(mtol)
    mtol = 1.0e-12;
end

% Compute element matrix at reference line
if (nargin < 6) || isempty(EM)
    EM = computeElementMatrix_refLine(trial_der, test_der, bs, mtol);
    if inverse
        EM = inv(EM);
    end
end

nb = size(EM, 1);
[re, ce, ve] = find(EM);
if isempty(ve)
    A = sparse([], [], [], msh.nLElems * nb, msh.nLElems * nb);
    return
end

% Assemble
h = msh.elemLength(:, msh.LElems);
J = msh.elemJac(:, msh.LElems);

J = J .* (2 ./ h).^(trial_der + test_der); 
if inverse
    J = 1 ./ J;
end
r = re + (0 : msh.nLElems - 1) * nb;
c = ce + (0 : msh.nLElems - 1) * nb;
v = J .* ve;
A = sparse(r, c, v, msh.nLElems * nb, msh.nLElems * nb);

end


