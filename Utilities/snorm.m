% Compute scaled discrete L2-norm, which is the discrete L2-norm divided by
% the square root of the number of elements in the vector.
% x : Vector for which scaled discrete L2-norm must be calculated
% Copyright (C) J. J. W. Van der Vget.

function res = snorm(v)

res = norm(v) / sqrt(length(v));
    
end