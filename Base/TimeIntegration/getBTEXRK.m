% set the Butcher tableau for explicit Runge-Kutta method
% type : 1 for the SSP RK method, 2 for the classical RK method, 
% Copyright (C) Xiangyi Meng.

function bt = getBTEXRK(order, type)

if (nargin < 1)
    error('Not enough arguments')
end

if (nargin < 2) || isempty(type)
    type = 2;
end

bt.type = type;
switch type
    % SSP RK method
    case 1
        switch order
            % Euler forward
            case 1
                bt.nstages = 1;
                bt.order = 1;
                bt.Alpha = 1;
                bt.Beta  = 1;
                bt.Gamma = 0;
                
            % second order
            case 2
                bt.nstages = 2;
                bt.order = 2;
                bt.Alpha = [1, 0; 1 / 2, 1 / 2];
                bt.Beta  = [1, 0;     0, 1 / 2];
                bt.Gamma = [0; 1];
                
            % third order
            case 3
                bt.nstages = 3;
                bt.order = 3;
                bt.Alpha = [1, 0, 0; 3 / 4, 1 / 4, 0; 1 / 3, 0, 2 / 3];
                bt.Beta  = [1, 0, 0;     0, 1 / 4, 0;     0, 0, 2 / 3];
                bt.Gamma = [0; 1; 1 / 2];
                
            % fourth order
            case 4
                bt.nstages = 5;
                bt.order = 4;
                bt.Alpha = [               1, 0, 0, 0, 0;0.44437049406734, 0.55562950593266, 0, 0, 0; 0.62010185138540, 0, 0.37989814861460, 0, 0; 0.17807995410773, 0, 0, 0.82192004589227, 0; 0.00683325884039, 0, 0.51723167208978, 0.12759831133288, 0.34833675773694];
                bt.Beta  = [0.39175222700392, 0, 0, 0, 0;               0, 0.36841059262959, 0, 0, 0;                0, 0, 0.25189177424738, 0, 0;                0, 0, 0, 0.54497475021237, 0;                0, 0,                0, 0.08460416338212, 0.22600748319395];
                bt.Gamma = [0; 0.39175222700392; 0.58607968896780; 0.47454236302687; 0.93501063100924];
            
            otherwise
                error('Not implemented SSP RK method')
        end
        
    % classical RK method
    case 2
        switch order
            % Euler forward
            case 1
                bt.nstages = 1;
                bt.order = 1;
                bt.A = 0;
                bt.b = 1;
                bt.c = 0;

            % second order
            case 2
                bt.nstages = 2;
                bt.order = 2;
                bt.A = [0, 0; 1, 0];
                bt.b = [0.5; 0.5];
                bt.c = [0; 1];

            % third order
            case 3
                bt.nstages = 3;
                bt.order = 3;
                bt.A = [0, 0, 0; 1, 0, 0; 0.25, 0.25, 0];
                bt.b = [1 / 6; 1 / 6; 2 / 3];
                bt.c = [0; 1; 0.5];

            % fourth order
            case 4
                bt.nstages = 4;
                bt.order = 4;
                bt.A = [  0,   0,   0, 0;
                        0.5,   0,   0, 0;
                          0, 0.5,   0, 0;
                          0,   0,   1, 0];
                bt.b = [1 / 6; 1 / 3; 1 / 3; 1 / 6];
                bt.c = [0; 0.5; 0.5; 1];

            % fifth order
            case 5
                bt.nstages = 6;
                bt.order = 5;
                bt.A = [       0,       0,         0,        0,          0, 0;
                             0.5,       0,         0,        0,          0, 0;
                            0.25,    0.25,         0,        0,          0, 0; 
                               0,      -1,         2,        0,          0, 0;
                          7 / 27, 10 / 27,         0,   1 / 27,          0, 0;
                        28 / 625,  -1 / 5, 546 / 625, 54 / 625, -378 / 625, 0];
                bt.b = [1 / 24; 0; 0; 5 / 48; 27 / 56; 125 / 336];
                bt.c = [0; 0.5; 0.5; 1; 2 / 3; 0.2];

            otherwise
                error('Not implemented classical RK method')
        end
           
    otherwise
        error('Not implemented type of the EXRK method')
end

