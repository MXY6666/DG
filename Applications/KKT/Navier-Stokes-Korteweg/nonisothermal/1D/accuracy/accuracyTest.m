%**************************************************************************
%                       problem description
%**************************************************************************
% Use DIRK-DG method to compute
% rho_t + (rho * u)_x = 0,
% (rho * u)_t + (rho * u^2 + p)_x = (tau + xi)_x,
% E_t + ((E + p) * u)_x = ((tau + xi) * u - q - jE)_x,
% where 
% E   = rho * (e + 1 / 2 * u^2) + 1 / 2 * 1 / We * rho_x^2,
% e   = 8 / 3 * Cv * T - 3 * rho,
% p   = 8 * rho * T / (3 - rho) - 3 * rho^2,
% tau = 1 / Re * 4 / 3 * u_x,
% xi  = 1 / We * (rho * rho_xx - 1 / 2 * rho_x^2),
% q   = -(8 * Cv) / (3 * We * Pr) * T_x,
% jE  = 1 / We * rho * rho_x * u_x.
% periodic boundary condition
% N            : number of elements
% tp           : time to print
% basisType    : type of basis functions
% DIRKType     : type of the DIRK time integration method
% inviscidFlux : 1 for Lax-Friedrichs flux
%                2 for local Lax-Friedrichs flux
% KKT          : parameters for the KKT-limiter
% newton       : parameters fot newton's method
% pat          : pat = 0 means only plot at the final time
%              : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = accuracyTest(N, tp, basisType, DIRKType, inviscidFlux, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-14;

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(basisType)
    basisType = 101;
end
if all(basisType ~= [100 : 104, 110 : 114, 120 : 124])
    error('Wrong basis type')
end

if (nargin < 4) || isempty(DIRKType)
    DIRKType = 11;
end 

if (nargin < 5) || isempty(inviscidFlux)
    inviscidFlux = 2;
end 

if (nargin < 8) || isempty(pat)
    pat = 0;
end 

% Choose test case
tc = createTestCase(50, 1000, 0.843, 5.375);

% quadrature rule
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 2, 102);

% basis function set data
bs = setBasisFunctionSet_line(quad1, basisType);

% mesh 
msh = setLineMesh_line(tc.dm, N, [1, 1], 101, 0, []);

% Augment mesh data
md = computeMeshData_line(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad1, bs);

% global matrix 
[IME, A, M, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i] = assembleMatrix(msh, md, tc, bs, mtol);

% initial numerical solution
[U0, W0, Lamb] = computeInitialSolution(msh, A, tc, quad1, bs, KKT, newton, mtol);

% CFL number
[CFL, CFLMax] = setCFLNumber;

% Butcher tableau
bt = getBTDIRK(k + 1, DIRKType);        

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial density
if (pat == 1)
    plotSols(msh, U0, W0, Lamb, 0, tc, quad1, bs, KKT)
end

% integrate in time using IMEX or SDC methods
t = 0;
nt = 0;
dt = setdt(msh, md, U0, W0, t, CFL, tp, tc, bs); 
while (t < tp - 1.0e-12)
    nt = nt + 1;       
    [U0, W0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, A, pr12, pc12, pv12, pr21, pc21, pv21, pr22i, pc22i, pv22i, IME, tc, quad1, bs, bt, inviscidFlux, KKT, newton, mtol); 
    
    % Plot the initial density
    if (pat == 1)
        plotSols(msh, U0, W0, Lamb, t, tc, quad1, bs, KKT)
    end
end 
RHO = reshape(U0(1 : end / 3), [bs.nb, msh.nLElems]);
M   = reshape(U0(end / 3 + 1 : 2 * end / 3), [bs.nb, msh.nLElems]);
E   = reshape(U0(2 * end / 3 + 1 : end), [bs.nb, msh.nLElems]);     
W = reshape(W0, [bs.nb, msh.nLElems]); 

% Plot 
if (pat == 0)
    plotSols(msh, U0, W0, Lamb, t, tc, quad1, bs, KKT)
end

% Compute error
quad1 = GaussQuadratureRule_line(7, quad1.type);
bs = setBasisFunctionSet_line(quad1, basisType);
[L2err, Lierr] = computeError(msh, RHO, M, E, W, t, tc, quad1, bs);

end





















        
        
        