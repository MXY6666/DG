% L-inf, L-2 errors and convergences of order
clc; clear
format long
%******************parameters setting**************************************
% problem index
pro = 3;
% number of elements
N = 200;  
% time to print
tp = 0.012;             
% type of basis functions
basisType = 102;      
% type of DIRK method
DIRKType = 11;
% slope limiter type
slo = 2;
% parameters for KKT-limiter
KKT.enable       = true;
KKT.conservation = false; % not implemented, keep false
KKT.positivity   = true;
KKT              = setKKTMin(KKT, 1.e-10);
KKT.Lamb0        = 1.e-6;
% parameters for newton method
newton.alpha                    = 1.e-4;
newton.alpha_bar                = 1.e-12;
newton.beta                     = 0.5;
newton.gamma                    = 0.5;
newton.sigma                    = 1.e-14;
newton.delta                    = 1.e-12;
newton.restart                  = 10;
newton.maxit                    = 50;
newton.maxitLineSearch          = 15;
newton.tolres                   = 1.e-8;
newton.toldV                    = 1.e-6;
newton.toldVLocal               = 1.e-13;
newton.nref                     = 4;
newton.scaleMat                 = false;
newton.linearSolver.type        = 1;
newton.linearSolver.refit       = 0;
newton.linearSolver.tolcor      = 1.e-15;
newton.linearSolver.precision   = 1;
newton.linearSolver.restart     = 5;
newton.linearSolver.maxit       = 100;
newton.linearSolver.tolres      = 1.e-10;
newton.linearSolver.ilu.type    = 'nofill'; % nofill, crout or ilutp
newton.linearSolver.ilu.droptol = 1.e-6;
newton.linearSolver.ilu.milu    = 'off'; % row, col or off
% Plot all the time or not
pat = 0;
%**************************************************************************
Riemann(pro, N, tp, basisType, DIRKType, slo, KKT, newton, pat) 











        
        
        
        