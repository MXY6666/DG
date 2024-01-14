%**************************************************************************
%                       problem description
%**************************************************************************
% Use DIRK-DG method to compute
% rho_t + (rho * u)_x + (rho * v)_y = 0,
% (rho * u)_t + (rho * u^2 + p)_x + (rho * u * v)_y = 0,
% (rho * v)_t + (rho * u * v)_x + (rho * v^2 + p)_y = 0,
% E_t + (u * (E + p))_x + (v * (E + p))_y = 0,
% where p = (gamma - 1) * (E - 1 / 2 * rho * (u^2 + v^2)).
% N         : number of elements in x and y direction
% tp        : time to print
% basisType : type of basis functions
% DIRKType  : type of the DIRK time integration method
% slo       : slo = 0 means no slope limiter
%             slo = 1 means TVD limiter 
%             slo = 2 means TVB limiter
% KKT       : parameters for the KKT-limiter
% newton    : parameters fot newton's method
% pat       : pat = 0 means only plot at the final time
%           : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function Sedov(N, tp, basisType, DIRKType, slo, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-14;

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(basisType)
    basisType = 201;
end
if (mod(floor(basisType / 10), 10) ~= 0)
    error('Only orthogonal basis functions are supported')
end

if (nargin < 4) || isempty(DIRKType)
    DIRKType = 11;
end 

if (nargin < 5) || isempty(slo)
    slo = 2;
end
if (slo ~= 0)
    if (mod(floor(basisType / 10), 10) ~= 0)
        error('Only orthogonal basis functions are supported when slope limiter is used')
    end
    
    if (mod(DIRKType, 10) ~= 1)
        error('Only standard formulation of the DIRK method is supported when slope limiter is used')
    end
end

if (nargin < 8) || isempty(pat)
    pat = 0;
end  

% Choose test case
tc = createTestCase(1.4);

% quadrature rule
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 1, 101);
quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);

% basis function set data
bs = setBasisFunctionSet_square(quad1, quad2, basisType);

% mesh
msh = setRectMesh_rect(tc.dm, N, [8, 7, 8, 7], 201, 0, [], []);

% Augment mesh data
md = computeMeshData_rect(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad2, bs);

% global mass matrix and inverse of mass matrix at reference square
M = assembleMatrix(msh, bs, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, md, tc, quad2, bs, slo, KKT, newton, mtol);

% CFL number
[CFL, CFLMax] = setCFLNumber;     

% Butcher tableau
bt = getBTDIRK(k + 1, DIRKType);

% Plot the initial density
if (pat == 1)
    plotSols(msh, U0, Lamb, 0, tc, quad2, bs, KKT)
end

%**************************************************************************
%                                 solve
%**************************************************************************
% integration in time using DIRK methods
t  = 0;
nt = 0;
dt = setdt(msh, md, U0, t, CFL, tp, tc, bs); 
while(t < tp - 1.0e-12)
    nt = nt + 1;
    [U0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, tc, quad2, bs, bt, slo, KKT, newton, mtol);      
    
    % Plot the initial density
    if (pat == 1)
        plotSols(msh, U0, Lamb, t, tc, quad2, bs, KKT)
    end 
end    

% Plot 
if (pat == 0)
    plotSols(msh, U0, Lamb, t, tc, quad2, bs, KKT)
end

end





















        
        
        