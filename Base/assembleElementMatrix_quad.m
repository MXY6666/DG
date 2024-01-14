% msh       : linear quadrilateral mesh in 2D
% trial_der : derivative order of trial basis functions
% test_der  : derivative order of test basis functions
% bs        : basis function set data
% inverse   : inverse the element matrix or not
% mtol      : tolerance for the matrix entries
% Copyright (C) Xiangyi Meng.

function A = assembleElementMatrix_quad(msh, trial_der, test_der, bs, inverse, mtol)

if (nargin < 5)
    error('Not enough arguments')
end    

if (msh.type ~= 211) && (msh.type ~= 212)
    error('Wrong mesh type')
end

if (nargin < 5) || isempty(inverse)
    inverse = false;
end

if (nargin < 6) || isempty(mtol)
    mtol = 1.0e-12;
end

if (bs.type == 200)
    if ~all(trial_der(:) == [0; 0]) || ~all(test_der(:) == [0; 0])
        A = sparse([], [], [], msh.nLElems * bs.nb, msh.nLElems * bs.nb);
        return; 
    end
end

% Compute element matrix
if all(trial_der(:) == [0; 0])
    if all(test_der(:) == [0; 0])
        EM = bs.phitw{1} * (repmat(bs.phi{1}, [1, msh.nLElems]) .* repelem(msh.elemJac(:, msh.LElems), 1, bs.nb));
    elseif all(test_der(:) == [1; 0])
        EM = bs.phitw{2} * (repmat(bs.phi{1}, [1, msh.nLElems]) .* repelem(msh.elemJxix(:, msh.LElems), 1, bs.nb)) ...
           + bs.phitw{3} * (repmat(bs.phi{1}, [1, msh.nLElems]) .* repelem(msh.elemJetax(:, msh.LElems), 1, bs.nb));
    elseif all(test_der(:) == [0; 1])
        EM = bs.phitw{2} * (repmat(bs.phi{1}, [1, msh.nLElems]) .* repelem(msh.elemJxiy(:, msh.LElems), 1, bs.nb)) ...
           + bs.phitw{3} * (repmat(bs.phi{1}, [1, msh.nLElems]) .* repelem(msh.elemJetay(:, msh.LElems), 1, bs.nb));
    else
        error('Not implemented derivative order of trial basis funtion')
    end    
elseif all(trial_der(:) == [1; 0])
    error('Not implemented case yet')
elseif all(trial_der(:) == [0; 1])
    error('Not implemented case yet')
else
    error('Not implemented derivative order of trial basis funtion')
end

if inverse
    EM = mychop(EM, mtol);
    EM = cell2mat(cellfun(@inv, mat2cell(EM, bs.nb, ones(1, msh.nLElems) * bs.nb), 'UniformOutput', false));
end
EM = mychop(EM, mtol);

% Assemble
re = repmat((1 : bs.nb)', [bs.nb, 1]);
ce = repelem((1 : bs.nb)', bs.nb);
A  = sparse(re + (0 : msh.nLElems - 1) * bs.nb, ce + (0 : msh.nLElems - 1) * bs.nb, EM, msh.nLElems * bs.nb, msh.nLElems * bs.nb);

end


