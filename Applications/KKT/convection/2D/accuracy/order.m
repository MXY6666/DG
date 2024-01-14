% L^2 and L^inf errors and convergences of order
clc; clear
format long
%******************parameters setting**************************************
% advection coefficient
beta = [1, 1];
% number of elements in x direction
Nx = [8 16 32 64]; 
% number of elements in y direction
Ny = [8 16 32 64]; 
% time to print
tp = 0.1; 
% type of basis functions
basisType = 201; 
% type of DIRK method
DIRKType = 11;
% parameters for KKT-limiter
KKT.enable       = true;
KKT.conservation = false; % not implemented, keep false
KKT.positivity   = false;
KKT.maximum      = true;
KKT.TVD          = false; % not implemented, keep false
KKT              = setKKTMinMax(KKT, 1.e-10, 1);
KKT.Lamb0        = 1.e-6;
% parameters for the Newton methods
newton.alpha                  = 1.e-4;
newton.alpha_bar              = 1.e-12;
newton.beta                   = 0.5;
newton.gamma                  = 0.5;
newton.sigma                  = 1.e-14;
newton.delta                  = 1.e-12; 
newton.restart                = 10;
newton.maxit                  = 35;
newton.maxitLineSearch        = 15;
newton.tolres                 = 1.e-10;
newton.toldV                  = 1.e-5;
newton.toldVLocal             = 1.e-13;
newton.nref                   = 4;
newton.scaleMat               = true;
newton.linearSolver.type      = 1;
newton.linearSolver.refit     = 0;
newton.linearSolver.tolcor    = 1.e-15;
newton.linearSolver.precision = 1;
newton.linearSolver.restart   = 5;
newton.linearSolver.maxit     = 100;
newton.linearSolver.tolres    = 1.e-7;
% Plot all the time or not
pat = 0;
%**************************************************************************
l = length(Nx);
Lierr = zeros(1, l); Lior = zeros(1, l - 1);
L2err = zeros(1, l); L2or = zeros(1, l - 1);
for s = 1 : l    
    [L2err(s), Lierr(s)] =  accuracyTest(beta, [Nx(s), Ny(s)], tp, basisType, DIRKType, KKT, newton, pat);    
end

for s = 1 : l - 1   
    L2or(s) = log(L2err(s) / L2err(s + 1)) / log(2);
    Lior(s) = log(Lierr(s) / Lierr(s + 1)) / log(2);        
end

L2err
Lierr
 
L2or
Lior



        
        
        
        