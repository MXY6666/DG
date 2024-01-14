%L-inf,L-2 errors and convergences of order
clear; clc
format long
%******************parameters setting**************************************
% number of elements
N = 100; 
% time to print
tp = 0.316; 
% type of basis functions
basisType = 103; 
% type of DIRK method
DIRKType = 11;
% parameters for KKT-limiter
KKT.enable       = true;
KKT.conservation = false; % not implemented, keep false
KKT.positivity   = true;
KKT.maximum      = true;
KKT.TVD          = false; % not implemented, keep false
KKT              = setKKTMinMax(KKT, 1.e-10, 1);
KKT.Lamb0        = 1.e-6;
% parameters for the Newton methods
newton.alpha                 = 1.e-4;
newton.alpha_bar             = 1.e-12;
newton.beta                  = 0.5;
newton.gamma                 = 0.5;
newton.sigma                 = 1.e-14;
newton.delta                 = 1.e-12; 
newton.restart               = 10;
newton.maxit                 = 60;
newton.maxitLineSearch       = 10;
newton.tolres                = 1.e-8;
newton.toldV                 = 1.e-6;
newton.toldVLocal            = 1.e-13;
newton.nref                  = 4;
newton.scaleMat              = false;
newton.linearSolver.type     = 1;
newton.linearSolver.refit    = 0;
newton.inearSolver.tolcor    = 1.e-15;
newton.inearSolver.precision = 1;
newton.linearSolver.restart  = 5;
newton.linearSolver.maxit    = 100;
newton.linearSolver.tolres   = 1.e-7;
% Plot all the time or not
pat = 1;
%**************************************************************************
travelingWave(N, tp, basisType, DIRKType, KKT, newton, pat);    


