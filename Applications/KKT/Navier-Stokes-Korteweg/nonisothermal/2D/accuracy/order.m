% L^2 and L^inf errors and convergences of order  
clc; clear
format long
%******************parameters setting**************************************
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
% 1 for LF, 2 for LLF
inviscidFlux = 1;
% parameters for KKT-limiter
KKT.enable       = true;
KKT.conservation = false;
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
newton.scaleMat                 = true;
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

l = length(Nx);
Lierr = zeros(3, l); Lior = zeros(3, l - 1);
L2err = zeros(3, l); L2or = zeros(3, l - 1);
for s = 1 : l   
    [L2err(:, s), Lierr(:, s)] = accuracyTest([Nx(s), Ny(s)], tp, basisType, DIRKType, inviscidFlux, KKT, newton, pat);   
end

for s = 1 : l - 1   
    L2or(:, s) = log(L2err(:, s) ./ L2err(:, s + 1)) / log(2);
    Lior(:, s) = log(Lierr(:, s) ./ Lierr(:, s + 1)) / log(2);        
end

fprintf('\n\n')
fprintf('L2 errors of density:\n\n')
fprintf('%16.8e', L2err(1, :))
fprintf('\n\nL2 orders of density:\n\n')
fprintf('%16.11f', L2or(1, :))
fprintf('\n\nL-inf errors of density:\n\n')
fprintf('%16.8e', Lierr(1, :))
fprintf('\n\nL-inf orders of density:\n\n')
fprintf('%16.11f', Lior(1, :))

fprintf('\n\nL2 errors of velocity:\n\n')
fprintf('%16.8e', L2err(2, :))
fprintf('\n\nL2 orders of velocity:\n\n')
fprintf('%16.11f', L2or(2, :))
fprintf('\n\nL-inf errors of velocity:\n\n')
fprintf('%16.8e', Lierr(2, :))
fprintf('\n\nL-inf orders of velocity:\n\n')
fprintf('%16.11f', Lior(2, :))

fprintf('\n\nL2 errors of pressure:\n\n')
fprintf('%16.8e', L2err(3, :))
fprintf('\n\nL2 orders of pressure:\n\n')
fprintf('%16.11f', L2or(3, :))
fprintf('\n\nL-inf errors of pressure:\n\n')
fprintf('%16.8e', Lierr(3, :))
fprintf('\n\nL-inf orders of pressure:\n\n')
fprintf('%16.11f', Lior(3, :))
fprintf('\n')



        
        
        
        