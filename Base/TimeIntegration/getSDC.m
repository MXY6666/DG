% SDC is to sotre the integral of k + 1 Lagrange interpolation basis functions
% in k intervals [tn,0(tn) tn,1], [tn,1 tn,2], [tn,2 tn,3]...
% [tn,k-1, tn,k(tn + 1)](ignoring the time step dt) 
% Copyright (C) Xiangyi Meng.

function SDC = getSDC(order)

SDC.order = order;
switch order 
        
    case {1, 2} 
        
        SDC.points  = [-1, 1];
        SDC.weights = [1 / 2, 1 / 2];
        
    case 3
        
        SDC.points  = [-1, 0, 1];
        SDC.weights = [5 / 24, 1 / 3, -1 / 24; ...
                      -1 / 24, 1 / 3, 5 / 24];
        
    case 4       
        
        SDC.points  = [-1, -sqrt(1 / 5), sqrt(1 / 5), 1];
        SDC.weights = [(11 + sqrt(5)) / 120,       (25 - sqrt(5)) / 120, (25 - 13 * sqrt(5)) / 120,  (sqrt(5) - 1) / 120; ...
                              -sqrt(5) / 60,           7 * sqrt(5) / 60,          7 * sqrt(5) / 60,        -sqrt(5) / 60; ...
                        (sqrt(5) - 1) / 120,  (25 - 13 * sqrt(5)) / 120,      (25 - sqrt(5)) / 120, (11 + sqrt(5)) / 120];
            
    case 5
        SDC.points  = [-1, -sqrt(3 / 7), 0, sqrt(3 / 7), 1];
        SDC.weights = [      (119 + 3 * sqrt(21)) / 1960,  (343 - 9 * sqrt(21)) / 2520, 8 / 45 - 32 / 35 / sqrt(21), (343 - 69 * sqrt(21)) / 2520,         3 * (sqrt(21) - 7) / 1960; ...
                       -3 * (105 + 8 * sqrt(21)) / 15680,         269 / 320 / sqrt(21),          32 / 35 / sqrt(21),         -61 / 320 / sqrt(21),   9 / 448 - 3 / 280 * sqrt(3 / 7); ...
                         9 / 448 - 3 / 280 * sqrt(3 / 7),         -61 / 320 / sqrt(21),          32 / 35 / sqrt(21),         269 / 320 / sqrt(21), -3 * (105 + 8 * sqrt(21)) / 15680; ...
                               3 * (sqrt(21) - 7) / 1960, (343 - 69 * sqrt(21)) / 2520, 8 / 45 - 32 / 35 / sqrt(21),  (343 - 9 * sqrt(21)) / 2520,       (119 + 3 * sqrt(21)) / 1960];
        
    otherwise
        error('Not implemented SDC time integration')
end

end
