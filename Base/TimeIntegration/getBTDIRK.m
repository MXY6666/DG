% set the Butcher tableau for implicit Runge-Kutta method
% type : 11 for strongly S-stable DIRK method in standard formulation 
%        12 for strongly S-stable DIRK method in Skvortsov formulation
%        21 for A-stable DIRK method in standard formulation 
%        22 for A-stable DIRK method in Skvortsov formulation
% Copyright (C) Xiangyi Meng.

function bt = getBTDIRK(order, type)

if (nargin < 1)
    error('Not enough arguments')
end

if (nargin < 2) || isempty(type)
    type = 22;
end

bt.type = type;
% formulation of the DIRK method, standard or Skvortsov
bt.form = mod(type, 10);

switch type
    % strongly S-stable DIRK method
    case {11, 12}
        switch order
            % Euler backward
            case 1
                bt.nstages = 1;
                bt.order = 1;
                bt.A = 1;
                bt.b = 1;
                bt.c = 1;
                bt.diagCoe = 1;
            
            % second order
            case 2
                alpha = 1 - 0.5 * sqrt(2);
                
                bt.nstages = 2;
                bt.order = 2;
                bt.A = [    alpha,     0;
                        1 - alpha, alpha];
                bt.b = [1 - alpha; alpha];
                bt.c = [alpha; 1];
                bt.diagCoe = alpha;
                
            % third order
            case 3
                gamma = 0.435866521508;
                b2 = 0.25 * (5 - 20 * gamma + 6 * gamma^2);
                
                bt.nstages = 3;
                bt.order = 3;
                bt.A = [          gamma,     0, 0;
                        (1 - gamma) / 2, gamma, 0;
                         1 - b2 - gamma,    b2, gamma];
                bt.b = [1 - b2 - gamma; b2; gamma];
                bt.c = [gamma; (1 + gamma) / 2; 1];      
                bt.diagCoe = gamma;
                
            % fourth order
            case 4
                bt.nstages = 5;
                bt.order = 4;
                bt.A = [1 / 4,     0,     0,       0,     0; 
                       -1 / 4, 1 / 4,     0,       0,     0;
                        1 / 8, 1 / 8, 1 / 4,       0,     0;
                       -3 / 2, 3 / 4, 3 / 2,   1 / 4,     0;
                            0, 1 / 6, 2 / 3, -1 / 12, 1 / 4];
                bt.b = [0; 1 / 6; 2 / 3; -1 / 12; 1 / 4];
                bt.c = [1 / 4; 0; 1 / 2; 1; 1];  
                bt.diagCoe = 1 / 4;
                
            % fifth order
            case 5
                tmp1 = sqrt(6);
                tmp2 = (6 - tmp1) / 10;
                
                bt.nstages = 5;    
                bt.order = 5;            
                bt.A = [                          tmp2,                              0,                          0,                        0,    0;
                                  (-6 + 5 * tmp1) / 14,                           tmp2,                          0,                        0,    0;
                             (888 + 607 * tmp1) / 2850,      (126 - 161 * tmp1) / 1425,                       tmp2,                        0,    0;
                          (3153 - 3082 * tmp1) / 14250,   (3213 + 1148 * tmp1) / 28500,   (-267 + 88 * tmp1) / 500,                     tmp2,    0;
                       (-32583 + 14638 * tmp1) / 71250, (-17199 + 364 * tmp1) / 142500, (1329 - 544 * tmp1) / 2500, (-96 + 131 * tmp1) / 625, tmp2];
                bt.b = [0; 0; 1 / 9; (16 - tmp1) / 36; (16 + tmp1) / 36];
                bt.c = [tmp2; (6 + 9 * tmp1) / 35; 1; (4 - tmp1) / 10; (4 + tmp1) / 10];                
                bt.diagCoe = tmp2;
            otherwise
                error('Not implemented strongly S-stable DIRK method')         
        end
        
    % A-stable DIRK method
    case {21, 22}
        switch order
            % Euler backward
            case 1
                bt.nstages = 1;
                bt.order = 1;
                bt.A = 1;
                bt.b = 1;
                bt.c = 1;
                bt.diagCoe = 1;
                
            % second order
            case 2
                bt.nstages = 1;
                bt.order = 2;
                bt.A = 0.5;
                bt.b = 1;
                bt.c = 0.5;
                bt.diagCoe = 0.5;

            % third order
            case 3
                gamma = 1 / 2 + sqrt(3) / 6 ;

                bt.nstages = 2;
                bt.order = 3;
                bt.A = [gamma, 0; 1 - 2 * gamma, gamma];
                bt.b = [0.5; 0.5];
                bt.c = [gamma, 1 - gamma];
                bt.diagCoe = gamma;

            % fourth order
            case 4
                gamma = 0.5 + cos(pi / 18) / sqrt(3);  
                delta = 1 / (6 * (2 * gamma - 1)^2);

                bt.nstages = 3;
                bt.order = 4;  
                bt.A = [      gamma,             0,     0; 
                        0.5 - gamma,         gamma,     0; 
                          2 * gamma, 1 - 4 * gamma, gamma];
                bt.b = [delta; 1 - 2 * delta; delta];
                bt.c = [gamma; 0.5; 1 - gamma];
                bt.diagCoe = gamma;

            otherwise
                error('Not implemented A-stable DIRK method')
        end 
    
    otherwise
        error('Not implemented type of the DIRK method')
    
end

switch bt.form
    case 1
        bt.Knstages = 1;
        bt.Lnstages = bt.nstages;
    case 2
        bt.Knstages = bt.nstages;
        bt.Lnstages = 1;
end

end

