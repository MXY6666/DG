%**************************************************************************
%                       problem description
%**************************************************************************
% Use EXRK-DG method to compute u_t - nu * u_xx + f(u) = 0
% N         : number of elements
% tp        : time to print
% basisType : type of basis functions
% flux      : a struct containing parameters for the flux
% DIRKType  : type of the DIRK time integration method
% KKT       : parameters for the KKT-limiter
% newton    : parameters for the Newton methods
% pat       : pat = 0 means only plot at the final time
%           : pat = 1 means plot all the time
% Copyright (C) Xiangyi Meng.

function travelingWave(N, tp, basisType, flux, DIRKType, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-12;

if (nargin < 2)
    error('Not enough arguments')
end

if (nargin < 3) || isempty(basisType)
    basisType = 101;
end
if all(basisType ~= [100 : 104, 111 : 114])
    error('Wrong basis type')
end

if (nargin < 4) || isempty(DIRKType)
    DIRKType = 11;
end 

if (nargin < 7) || isempty(pat)
    pat = 0;
end  

% Create test case
tc = createTestCase;

% quadrature rule
k = mod(basisType, 10);
quad = GaussQuadratureRule_line(k + 2, 102);

% basis function set data
bs = setBasisFunctionSet_line(quad, basisType);

% mesh 
msh = setLineMesh_line(tc.dm, N, [2, 2], 101, 0, []);

% Augment mesh data
md = computeMeshData_line(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad, bs);

% global matrices
[IME, M, A1, S2, S] = assembleMatrix(msh, md, flux, tc, bs, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, tc, quad, bs, KKT, newton);

% CFL number
cfl = setCFLNumber(k);          

% Butcher tableau for DIRK time integration method
bt = getBTDIRK(k + 1, DIRKType);

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial solution
if (pat == 1)
    Q0 = computeAuxiliaryVariable(msh, md, U0, 0, A1, IME, tc, bs);
    plotSols(msh, U0, Q0, Lamb, 0, tc, quad, bs, KKT);
end

% integration in time using DIRK method
nt = 0;
t  = 0;
dt = setdt(msh, t, cfl, tp);
while (t < tp - 1.0e-12)
    nt = nt + 1;
    maxdt = setdt(msh, t, cfl, tp); 
    [U0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, maxdt, M, A1, S2, S, IME, tc, quad, bs, bt, KKT, newton, mtol);        
    
    % Plot
    if (pat == 1)
        Q0 = computeAuxiliaryVariable(msh, md, U0, 0, A1, IME, tc, bs);
        plotSols(msh, U0, Q0, Lamb, t, tc, quad, bs, KKT);
    end    
end    

% Plot
if (pat == 0)
    Q0 = computeAuxiliaryVariable(msh, md, U0, 0, A1, IME, tc, bs);
    plotSols(msh, U0, Q0, Lamb, t, tc, quad, bs, KKT);   
end 

end





















        
        
        