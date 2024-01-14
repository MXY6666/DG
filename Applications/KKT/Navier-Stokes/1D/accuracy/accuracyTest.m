%**************************************************************************
%                       problem description
%**************************************************************************
% Use DIRK-DG method to compute
% rho_t + (rho * u)_x = 0,
% (rho * u)_t + (rho * u^2 + p)_x = tau_x,
% E_t + (u * (E + p))_x = (tau * u - q)_x,
% where p = (gamma - 1) * rho * e, 
% E = rho * (e + 0.5 * u^2),
% e = cv * T,
% tau = mu * 4 / 3 * u_x,
% q = -k * T,
% periodic boundary condition
% N         : number of elements
% tp        : time to print
% basisType : type of basis functions
% DIRKType  : type of the DIRK time integration method
% KKT       : parameters for the KKT-limiter
% newton    : parameters fot newton's method
% pat       : pat = 0 means only plot at the final time
%           : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = accuracyTest(N, tp, basisType, DIRKType, KKT, newton, pat)        
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

if (nargin < 7) || isempty(pat)
    pat = 0;
end 

% Choose test case
tc = createTestCase(1.4, 0.7, 200);

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
[IME, M, pr12, pc12, pv12] = assembleMatrix(msh, md, tc, bs, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, tc, quad1, bs, KKT, newton, mtol);

% CFL number
[CFL, CFLMax] = setCFLNumber;

% Butcher tableau
bt = getBTDIRK(k + 1, DIRKType);

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial density
if (pat == 1)
    plotSols(msh, U0, Lamb, 0, tc, quad1, bs, KKT)
end

% integration in time using DIRK method
nt = 0;
t  = 0;
dt = setdt(msh, md, U0, t, CFL, tp, tc, bs); 
while (t < tp - 1.0e-12)
    nt = nt + 1;        
    [U0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, pr12, pc12, pv12, IME, tc, quad1, bs, bt, KKT, newton, mtol);
    
    % Plot the initial density
    if (pat == 1)
        plotSols(msh, U0, Lamb, t, tc, quad1, bs, KKT)
    end
end 
RHO = reshape(U0(1 : end / 3), [bs.nb, msh.nLElems]);
M   = reshape(U0(end / 3 + 1 : 2 * end / 3), [bs.nb, msh.nLElems]);
E   = reshape(U0(2 * end / 3 + 1 : end), [bs.nb, msh.nLElems]);     

% Plot 
if (pat == 0)
    plotSols(msh, U0, Lamb, t, tc, quad1, bs, KKT)
end

% Compute error
quad1 = GaussQuadratureRule_line(7, quad1.type);
bs = setBasisFunctionSet_line(quad1, basisType);
[L2err, Lierr] = computeError(msh, RHO, M, E, tp, tc, quad1, bs);

end





















        
        
        