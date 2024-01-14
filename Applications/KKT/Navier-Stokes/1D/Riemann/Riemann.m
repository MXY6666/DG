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
% pro : choose which problem to solve
%       pro = 1 (Sod problem, left rarefaction, middle contact discontinuity, right shock): 
%               rhoL =   1.0; uL = 0.0; pL = 1.0;
%               rhoR = 0.125; uR = 0.0; pR = 0.1;
%               tp = 0.25;
%       pro = 2 (Lax problem, left rarefaction, middle contact discontinuity, right shock) :
%               rhoL = 0.445; uL = 0.698; pL = 3.528;
%               rhoR = 0.5;   uR = 0;     pR = 0.571;
%               tp = 0.13;
%       pro = 3 (left rarefaction, middle contact discontinuity, right shock):
%               rhoL = 1.0; uL = 0.0; pL = 1000.0;
%               rhoR = 1.0; uR = 0.0; pR = 0.01; 
%               tp = 0.012;
%       pro = 4 (left shock, middle contact discontinuity, right rarefaction):
%               rhoL = 1.0; uL = 0.0; pL = 0.01;
%               rhoR = 1.0; uR = 0.0; pR = 100.0;
%               tp = 0.035;
%       pro = 5 (left shock, middle contact discontinuity, right shock):
%               rhoL = 5.99924; uL =  19.5975; pL = 460.894;
%               rhoR = 5.99242; uR = -6.19633; pR = 46.0950;
%               tp = 0.035;
%       pro = 6 (123 problem, two strong rarefactions and a stationary contact discontinuity):
%               rhoL = 1.0; uL = -2.0; pL = 0.4;
%               rhoR = 1.0; uR =  2.0; pR = 0.4; 
%               tp = 0.15;
%       pro = 7 (123 problem, two strong rarefactions and a stationary contact discontinuity):
%               rhoL = 7; uL = -1; pL = 0.2;
%               rhoR = 7; uR =  1; pR = 0.2; 
%               tp = 0.6;
%       pro = 8 (Leblanc problem):
%               rhoL = 10^4; uL = 0; pL = 10^4;
%               rhoR = 1;    uR = 0; pR = 1; 
%               tp = 0.12;
%       pro = 9 (Leblanc problem):
%               rhoL = 2;     uL = 0; pL = 10^9;
%               rhoR = 1.e-3; uR = 0; pR = 1; 
%               tp = 1.e-4;
%       pro = 10 (Sedov blast wave problem):
%               rho = 1; u = 0; 
%               p = 10^2 in middle element, 1.e-9 otherwise;
%               tp = 0.5;
%       pro = 11 (Sedov blast wave problem):
%               rho = 1; u = 0; 
%               E = 3.2e6 / h in middle element, 1.e-9 otherwise;
%               tp = 1.e-3;
% N         : number of elements
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

function Riemann(pro, N, tp, basisType, DIRKType, slo, KKT, newton, pat)        
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-14;

if (nargin < 3)
    error('Not enough arguments')
end

if all(pro ~= 1 : 11)
    error('Wrong problem index')
end
if (pro == 8)
    if (isinteger(0.3 * N) ~= 0)
        error('The number of elements multiplies 0.3 should be an integer for this problem')
    end
else
    if (mod(N, 2) ~= 0)
        error('The number of elements should be an even number for this problem')
    end
end

if (nargin < 4) || isempty(basisType)
    basisType = 101;
end
if (mod(floor(basisType / 10), 10) ~= 0)
    error('Only orthogonal basis functions are supported')
end

if (nargin < 5) || isempty(DIRKType)
    DIRKType = 11;
end 

if (nargin < 6) || isempty(slo)
    slo = 1;
end 
if (slo ~= 0 && mod(DIRKType, 10) ~= 1)
    error('Only standard formulation of the DIRK method is supported when slope limiter is used')
end

if (nargin < 9) || isempty(pat)
    pat = 0;
end 

% Choose test case
tc = createTestCase(pro, 1.4, 0.72, 2000);

% quadrature rule
k = mod(basisType, 10);
quad1 = GaussQuadratureRule_line(k + 2, 102);

% basis function set data
bs = setBasisFunctionSet_line(quad1, basisType);

% mesh 
msh = setLineMesh_line(tc.dm, N, [2, 2], 101, 0, []);

% Augment mesh data
md = computeMeshData_line(msh);

% Set DOFs for KKT-limiter
KKT = setKKTDOF(KKT, msh, quad1, bs);

% global matrix 
[IME, M, pr12, pc12, pv12] = assembleMatrix(msh, md, tc, bs, mtol);

% initial numerical solution
[U0, Lamb] = computeInitialSolution(msh, tc, quad1, bs, KKT, newton, mtol);

% CFL number
[CFL, CFLMax] = setCFLNumber(tc);

% Butcher tableau
bt = getBTDIRK(k + 1, DIRKType);

% points to plot
pn = setPlotPoints(msh, 1, [], quad1);

%**************************************************************************
%                                 solve
%**************************************************************************
% Plot the initial density
if (pat == 1)
    plotSols(msh, U0, Lamb, 0, pn, tc, quad1, bs, KKT)
end

% integration in time using DIRK method
nt = 0;
t  = 0;
dt = setdt(msh, md, U0, t, CFL, tp, tc, bs); 
while (t < tp - 1.0e-12)
    nt = nt + 1;       
    [U0, Lamb, t, dt] = computeOneTimeStepDIRK(msh, md, U0, Lamb, nt, t, dt, CFLMax, tp, M, pr12, pc12, pv12, IME, tc, quad1, bs, bt, slo, KKT, newton, mtol);
    
    % Plot the initial density
    if (pat == 1)
        plotSols(msh, U0, Lamb, t, pn, tc, quad1, bs, KKT)
    end
end 
   
% Plot 
if (pat == 0)
    plotSols(msh, U0, Lamb, t, pn, tc, quad1, bs, KKT)
end

end





















        
        
        