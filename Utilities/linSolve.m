% Solve linear system using iterative refinement given an LU decomposition
% of the matrix A computed with Matlab function lu
% A         : sparse matrix of linear system
% b         : righthand side of linear system 
% solver    : 1 for permuted LU-decomposition
%           : 2 for minimum-norm solution
%           : 3 for gmres solution
% Copyright (C) J. J. W. Van der Vget.

function x = linSolve(A, b, solver)

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(solver)
    solver.type      = 1;
    solver.refit     = 0;
    solver.tolcor    = 1.e-15;
    solver.precision = 1;
end
  
if (solver.type == 1 || solver.refit > 0)
    [L, U, P, Q, D] = lu(A);
end

% Compute solution of linear system
if (solver.type == 1)
    % Use permuted LU-decomposition (not always stable)
    % Turn off warnings for nearly singular matrix. They occur since the
    % Jacobian is sometimes rank deficit, but this does not affect the
    % Newton method.

    warning('off')
    x = Q * (U \ (L \ (P * (D \ b))));
    warning('on')

elseif (solver.type == 2)

    x = lsqminnorm(A, b);
    
elseif (solver.type == 3)
    
    if isfield(solver, 'ilu')
        [L, U] = ilu(A, solver.ilu);
    else
        L = []; U = [];
    end
    
    if isfield(solver, 'x0')
        x = gmres(A, b, solver.restart, solver.tolres, solver.maxit, L, U, solver.x0);
    else
        x = gmres(A, b, solver.restart, solver.tolres, solver.maxit, L, U);
    end

end

if (solver.refit > 0)
    % Perform iterative refinement of solution of linear system

    if (solver.precision == 1)
        % Perform iterative refinement using a residual compution in
        % double solver.precision. Note, this method is faster, but can only
        % give a rather limited improvement

        x = iterativeRefinementDoubleprecision(A, L, U, P, Q, D, x, b, solver.refit, solver.tolcor);

    elseif (solver.precision == 2)
        % Perform iterative refinement using a residual compution in
        % quadruple solver.precision. 

        x = iterativeRefinementQuadprecision(A, L, U, P, Q, D, x, b, solver.refit, solver.tolcor);

    end
end
    
end

