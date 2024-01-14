%**************************************************************************
%                       problem description
%**************************************************************************
% Use EXRK-LDG or SDC-LDG method to compute the isothermal Euler-Korteweg 
% equations with constant capillarity in 1D
% rho_t + (rho * u)_x = 0,
% (rho * u)_t + (rho * u^2)_x + rho * F'(rho)_x = kappa * rho * rho_xxx,
% periodic boundary condition
% N         : number of elements
% tp        : time to print
% basisType : type of basis functions
% flux      : a struct containing parameters for the flux
% KKT       : parameters for the KKT-limiter
% newton    : parameters for the Newton methods
% pat       : pat = 0 means only plot at the final time
%           : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function [L2err, Lierr] = accuracyTest(N, tp, basisType, flux, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-18;

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(basisType)
    basisType = 101;
end
if all(basisType ~= [100 : 104, 111 : 114])
    error('Wrong basis type')
end

if (nargin < 4) || isempty(flux)
    flux.alpha1 = 0;
    flux.alpha2 = 0;
end 

if (nargin < 7) || isempty(pat)
    pat = 0;
end  

% Create test case
tc = createTestCase(0.01, 0.85, 8/27, 1, 1);

% quadrature rule
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 1, 101);

% basis function set data
bs = setBasisFunctionSet_line(quad1, basisType);

% mesh 
msh = setLineMesh_line(tc.dm, N, [1, 1], 101, 0, []);

% Augment mesh data
md = computeMeshData_line(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad1, bs);

% inverse of mass matrix at reference line
ME  = computeElementMatrix_refLine(0, 0, bs);
IME = mychop(inv(ME), mtol);

% global matrix 
[AE, A2] = assembleMatrix(msh, md, flux, tc, bs, IME, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, tc, quad1, bs, KKT, newton);

% CFL number
cfl = setCFLNumber(k);

% time integration method
SDC = getSDC(k + 1);
SDC.theta = 1;

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial density
if (pat == 1)
    plotDensity(msh, U0, Lamb, 0, tc, quad1, bs, KKT)
end

% integration in time
t  = 0;   
nt = 0;
dt = setdt(msh, U0, t, tc, quad1, bs, cfl, tp, false);
while(t < tp - 1.0e-12)
    nt = nt + 1;
    maxdt = setdt(msh, U0, t, tc, quad1, bs, cfl, tp, false);
    [U0, Lamb, t, dt] = computeOneTimeStepSDC(msh, md, U0, Lamb, nt, t, dt, maxdt, flux, AE, A2, ME, IME, tc, quad1, bs, SDC, KKT, newton);  
    
    % Plot the initial density
    if (pat == 1)
        plotDensity(msh, U0, Lamb, t, tc, quad1, bs, KKT)
    end
end
U0 = reshape(U0, [bs.nb, 2 * msh.nLElems]);

% Plot 
if (pat == 0)
    plotDensity(msh, U0, Lamb, t, tc, quad1, bs, KKT)
end

% Compute error
quad1 = GaussQuadratureRule_line(10, 101);
bs = setBasisFunctionSet_line(quad1, basisType);
[L2err, Lierr] = computeError(msh, U0, tp, tc, quad1, bs);



end




