% Scale a matrix using row-column scaling matrix 
% Algorithm: A Symmetry Preserving Algorithm for Matrix Scaling, P.
% Knight, D. Ruiz, B. Ucar, INRIA Report 7552, 2012,  https://hal.inria.fr/inria-00569250v3
% Copyright (C) Xiangyi Meng.

function [A, PR, PC] = rowcolScale(A0, tol, maxit)

[m, n] = size(A0);
PR = speye(m);
PC = speye(n);
A  = A0;

for k = 1 : maxit
    % maximun norm of each row and col vector
    r = max(abs(A'));
    c = max(abs(A));
    r(r < eps) = 1;
    c(c < eps) = 1;
    
    % row-column scaling
    PR = PR * sparse(1 : m, 1 : m, 1 ./ sqrt(r), m, m);
    PC = PC * sparse(1 : n, 1 : n, 1 ./ sqrt(c), n, n);
    A  = PR * A0 * PC;
    
    if (max(abs(1 - r)) < tol && max(abs(1 - c)) < tol)
       break; 
    end
end

end