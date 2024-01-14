% Copyright (C) Xiangyi Meng.

function bt = getBTIMEX(order)

%some constants of the implicit-explicit Runge-Kutta Butcher tableau
switch order
    % stage = 1, order = 1
    case 1
        bt.AHat = [0, 0; 1, 0];
        bt.A = [0, 0; 0, 1];
        bt.bHat = [1; 0];
        bt.b = [0; 1];
        bt.c = [0; 1];
        bt.diagCoe = 1;
    
    % stage = 2， order = 2
    case 2
        gamma = 1 - sqrt(2) / 2;
        delta = 1 - 1 / (2 * gamma);
        
        bt.AHat = [    0,         0, 0;
                   gamma,         0, 0;
                   delta, 1 - delta, 0];
        bt.A = [0,         0,     0;
                0,     gamma,     0; 
                0, 1 - gamma, gamma];
        bt.bHat = [delta; 1 - delta; 0];
        bt.b = [0; 1 - gamma; gamma];
        bt.c = [0; gamma; 1];
        bt.diagCoe = gamma;
        
    % stage = 3， order = 3
    case 3
        gamma = 0.435866521508459;
        beta1 = -3 / 2 * gamma^2 + 4 * gamma - 1 / 4;
        beta2 = 3 / 2 * gamma^2 - 5 * gamma + 5 / 4;
        alpha1 = -0.35;
        alpha2 = (1 / 3 - 2 * gamma^2 - 2 * beta2 * alpha1 * gamma) / (gamma * (1 - gamma));
        
        bt.AHat = [                       0,          0,      0, 0;
                                      gamma,          0,      0, 0;
                   (1 + gamma) / 2 - alpha1,     alpha1,      0, 0;
                                          0, 1 - alpha2, alpha2, 0];
        bt.A = [0,               0,     0,     0;
                0,           gamma,     0,     0;
                0, (1 - gamma) / 2, gamma,     0;
                0,           beta1, beta2, gamma];
        bt.bHat = [0; beta1; beta2; gamma];
        bt.b = [0; beta1; beta2; gamma];
        bt.c = [0; gamma; (1 + gamma) / 2; 1];
        bt.diagCoe = gamma;
        
    otherwise
        error('Not implemented implicit-explicit Runge-Kutta time integration')
end

end

