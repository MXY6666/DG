%**************************************************************************
%                       problem description
%**************************************************************************
% Use EXRK_DG method to compute u_t + beta1 * u_x + beta1 * u_y = 0
% periodic boundary condition
% beta      : advection coefficient 
% N         : number of elements in x and y direction
% tp        : time to print
% basisType : type of basis functions
% DIRKType  : type of the DIRK time integration method
% KKT       : parameters for the KKT-limiter
% newton    : parameters for the Newton methods
% pat       : pat = 0 means only plot at the final time
%           : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function travelingWave(beta, N, tp, basisType, DIRKType, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-10;

if (nargin < 3)
    error('Not enough arguments')
end

if (nargin < 4) || isempty(basisType)
    basisType = 201;
end
if all(basisType ~= [200 : 204, 210 : 214, 220 : 224]) 
    error('Wrong basis type')
end

if (nargin < 5) || isempty(DIRKType)
    DIRKType = 11;
end

if (nargin < 8) || isempty(pat)
    pat = 0;
end  

% Choose test case
tc = createTestCase(beta, KKT);

% quadrature rule
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 1, 101);
quad2 = GaussQuadratureRule_square([k + 1, k + 1], 201);

% basis function set data
bs = setBasisFunctionSet_square(quad1, quad2, basisType);

% mesh 
msh = setRectMesh_rect(tc.dm, N, ones(1, 4), 201, 0, [], []);

% Augment mesh data
md = computeMeshData_rect(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad2, bs);

% global matrix
[M, S] = assembleMatrix(msh, md, tc, bs, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, tc, quad2, bs, KKT, newton);

% CFL number
cfl = setCFLNumber(k);   

% Butcher tableau for DIRK time integration method
bt = getBTDIRK(k + 1, DIRKType);

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial solution
if (pat == 1)
    plotSols(msh, U0, Lamb, 0, tc, quad2, bs, KKT);
end

% integration in time using DIRK method
nt    = 0;
t     = 0;
dt    = setdt(msh, t, tc, cfl, tp);
maxdt = dt;
while (t < tp - 1.0e-12)
    nt = nt + 1;
    [U0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, U0, Lamb, nt, t, dt, maxdt, M, S, tc, quad2, bs, bt, KKT, newton);        
    if (t + dt) > tp
        dt = tp - t;
    end  
    
    % Plot
    if (pat == 1)
        plotSols(msh, U0, Lamb, t, tc, quad2, bs, KKT);
    end     
end 

% Plot
if (pat == 0)
    plotSols(msh, U0, Lamb, t, tc, quad2, bs, KKT);   
end 

end



    