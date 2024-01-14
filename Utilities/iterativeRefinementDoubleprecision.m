% Improve the accuracy of the solution of the linear system Ax=b using iterative
% refinement. 
%
% Note, iterative refinement is useful to get a more accurate solution for 
% ill-conditioned linear systems.
%
% Note, iterative refinement requires that the residual is computed with 
% quadruple precision. In this function the residual is computed with only 
% double precision, which strongly limits the usefulness of iterative 
% refinement. A better alternative is provided by the function 
% iterativeRefinementQuadprecision  
%
% A                       Matrix of linear system
% L                       Lower triangular matrix in LU-decomposition
% U                       Upper triangular matrix in LU decomposition
% P                       Permutation matrix in sparse LU-decomposition
% Q                       Permutation matrix in sparse LU-decomposition
% D                       Scaling factor in sparse LU-decomposition
% x                       Solution of linear system Ax=b which needs to be
%                         improved using iterative refinement
% b                       Righthand vector of linear system
% 
% Copyright (C) J. J. W. Van der Vget.

function x = iterativeRefinementDoubleprecision(A, L, U, P, Q, D, x, b, refit, tolcor)

% Determine size of matrix A
[NA1, NA2] = size(A);

% Check if matrix A is square
if(NA1 ~= NA2)
    error('Iterative refinement requires a square matrix')
end

% Determine size of right hand side vector b
[Nb1, Nb2] = size(b);

% Check if size of vector b matches with matrix A
if(NA2 ~= Nb1)
    error('Vector b must match dimension of matrix A')
end

if(Nb2 ~= 1)
    error('Vector b in iterative refinement function can have only one column')
end

% Apply iterative refinement for matrices with a poor condition number
for n = 1 : refit

    % 1. Compute residual (should ideally be done in double precision)
    res = b - A * x;

    % Compute correction to solution x
    x_cor = Q * (U \ (L \ (P * (D \ res))));

    % Apply iterative refinement correction
    x = x + x_cor;

    % Compute maximum correction in iterative refinement
    x_cormax = max(abs(x_cor));

    if(x_cormax < tolcor)
        % Stop iterative refinement loop if correction is negligible
        break
    end

end
    
end