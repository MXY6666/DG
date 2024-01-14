% Improve the accuracy of the solution of the linear system Ax=b using iterative
% refinement. 
%
% Note, iterative refinement is useful to get a more accurate solution for 
% ill-conditioned linear systems.
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

function x = iterativeRefinementQuadprecision(A, L, U, P, Q, D, x, b, refit, tolcor)
 
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

% Find row and column indices of sparse matrix A
[iA, jA, vA] = find(A);

% Determine number of nonzero elements in each column of A
isA = zeros(NA2, 1);
ieA = zeros(NA2, 1);

isA(1) = 1;
ieA(1) = 1;
jc = 1;

for i = 2 : length(iA)
    if(jA(i) == jA(i - 1))
        ieA(jc) = i;
    else
        jc = jc + 1;
        isA(jc) = i;           
    end
end

ieA(NA2) = length(iA);

% Split A into an upper and lower part
cA = 134217729 * vA; 
hA = cA - (cA - vA);
lA = vA - hA;

% Apply iterative refinement for matrices with a poor condition number
for n = 1 : refit
    % Split x into an upper and lower part
    cx = 134217729 * x; 
    hx = cx - (cx - x);
    lx = x - hx;

    % 1. Compute residual (should ideally be done in double precision)
    % Initialize residual 
    res = zeros(NA1, 1);

    % Compute residual res=A*x-b
    err = zeros(NA1, 1);

    for j = 1 : NA2
        % Get non-zero elements sparse matrix A           
        i1 = isA(j);
        i2 = ieA(j);

        Ax = vA(i1 : i2) * x(j);

        lAtmp = lA(i1 : i2);
        hAtmp = hA(i1 : i2);
        errres = lAtmp * lx(j) - (((Ax - hAtmp * hx(j)) - lAtmp * hx(j)) - hAtmp * lx(j));          
           
        restmp = res(iA(i1 : i2));
        xs = restmp + Ax;
        zs = xs - restmp;
        errres = (restmp - (xs - zs)) + (Ax - zs) + errres;

        res(iA(i1 : i2)) = xs;  %Correct??
        err(iA(i1 : i2)) = err(iA(i1 : i2)) + errres;  %Correct??
    end

    % Subtract righthand side vector b
    [res, errsum] = twoSumV(res, -b);

    % Update error with contribution from summation
    err = err + errsum;

    % Correct residual for error in residual computation
    res = res + err;

    % Compute iterative refinement correction
    x_cor = Q * (U \ (L \ (P * (D \ (-res)))));

    % Apply iterative refinement correction to solution x
    x = x + x_cor;

    % Compute maximum size correction
    x_cor_max = max(abs(x_cor));

    if(x_cor_max < tolcor)
        % Stop iterative refinement loop if correction is negligible
        break
    end
end
    
end