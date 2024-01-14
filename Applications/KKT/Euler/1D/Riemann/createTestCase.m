function tc = createTestCase(pro, gamma)

% problem index
tc.pro = pro;

% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;

switch pro
    % Sod problem, left rarefaction, middle contact discontinuity, right shock
    case 1
        tc.rhoL = 1.0;   tc.uL = 0.0; tc.pL = 1.0;
        tc.rhoR = 0.125; tc.uR = 0.0; tc.pR = 0.1; 
        
    % Lax problem, left rarefaction, middle contact discontinuity, right shock
    case 2
        tc.rhoL = 0.445; tc.uL = 0.698; tc.pL = 3.528;
        tc.rhoR = 0.5;   tc.uR = 0;     tc.pR = 0.571;
                   
    % left rarefaction, middle contact discontinuity, right shock
    case 3
        tc.rhoL = 1.0; tc.uL = 0.0; tc.pL = 1000.0;
        tc.rhoR = 1.0; tc.uR = 0.0; tc.pR = 0.01;  
        
    % left shock, middle contact discontinuity, right rarefaction
    case 4
        tc.rhoL = 1.0; tc.uL = 0.0; tc.pL = 0.01;
        tc.rhoR = 1.0; tc.uR = 0.0; tc.pR = 100.0;   
       
    % left shock, middle contact discontinuity, right shock
    case 5
        tc.rhoL = 5.99924; tc.uL = 19.5975;  tc.pL = 460.894;
        tc.rhoR = 5.99242; tc.uR = -6.19633; tc.pR = 46.0950; 
         
    % 123 problem, two strong rarefactions and a stationary contact discontinuity
    case 6
        tc.rhoL = 1.0; tc.uL = -2.0; tc.pL = 0.4;
        tc.rhoR = 1.0; tc.uR = 2.0;  tc.pR = 0.4;
        
    % 123 problem, two strong rarefactions and a stationary contact discontinuity
    case 7
        tc.rhoL = 7; tc.uL = -1; tc.pL = 0.2;
        tc.rhoR = 7; tc.uR =  1; tc.pR = 0.2; 
        
    % Leblanc problem
    case 8
        tc.rhoL = 10^4; tc.uL = 0; tc.pL = 10^4;
        tc.rhoR = 1;    tc.uR = 0; tc.pR = 1;  
        
    % Leblanc problem
    case 9
        tc.rhoL = 2;     tc.uL = 0; tc.pL = 10^9;
        tc.rhoR = 1.e-3; tc.uR = 0; tc.pR = 1; 
        
    % Sedov blast wave problem
    case 10
        tc.rhoL = 1; tc.uL = 0; tc.pL = 1.e-9;
        tc.rhoR = 1; tc.uR = 0; tc.pR = 1.e-9;
        tc.pm = 10^2;
        
    % Sedov blast wave problem
    case 11 
        tc.rhoL = 1; tc.uL = 0; tc.EL = 1.e-9;
        tc.rhoR = 1; tc.uR = 0; tc.ER = 1.e-9;
        tc.Em = 3.2e6;
             
    otherwise
        error('Not implemented problem index')
end

tc.mL = tc.rhoL * tc.uL; 
tc.mR = tc.rhoR * tc.uR; 
if (pro == 10)
    tc.EL = tc.pL / (tc.gamma - 1);
    tc.ER = tc.pR / (tc.gamma - 1);
    tc.Em = tc.pm / (tc.gamma - 1);
elseif (pro == 11)
    tc.pL = (tc.gamma - 1) * tc.EL;
    tc.pR = (tc.gamma - 1) * tc.ER;
    tc.pm = (tc.gamma - 1) * tc.Em;
else   
    tc.EL = tc.pL / (tc.gamma - 1) + tc.rhoL * tc.uL^2 / 2;
    tc.ER = tc.pR / (tc.gamma - 1) + tc.rhoR * tc.uR^2 / 2;

    % initial solution
    tc.rho0 = @(x)tc.rhoL * (x < 0) + tc.rhoR * (x >= 0);
    tc.u0   = @(x)tc.uL   * (x < 0) + tc.uR   * (x >= 0);
    tc.p0   = @(x)tc.pL   * (x < 0) + tc.pR   * (x >= 0);
    tc.m0   = @(x)tc.mL   * (x < 0) + tc.mR   * (x >= 0);
    tc.E0   = @(x)tc.EL   * (x < 0) + tc.ER   * (x >= 0);

    % some constants required to compute p-star and exact solution
    tc.aL = sqrt(tc.gamma * tc.pL ./ tc.rhoL); tc.aR = sqrt(tc.gamma * tc.pR ./ tc.rhoR);
    tc.AL = 2 / ((tc.gamma + 1) * tc.rhoL); tc.BL = (tc.gamma - 1) / (tc.gamma + 1) * tc.pL;
    tc.AR = 2 / ((tc.gamma + 1) * tc.rhoR); tc.BR = (tc.gamma - 1) / (tc.gamma + 1) * tc.pR;

    %***************some possible choices of initial guess to p-star*****************************
    tc.TOL = 1e-6;

    % Two-Rarefaction approximation
    tc.pTR = ((tc.aL + tc.aR - (tc.gamma - 1) * (tc.uR - tc.uL) / 2) / (tc.aL / tc.pL^((tc.gamma - 1) / (2 * tc.gamma)) + tc.aR / tc.pR^((tc.gamma - 1) / (2 * tc.gamma))))^(2 * tc.gamma / (tc.gamma - 1));

    % linearised solution
    tc.pPV = (tc.pL + tc.pR) / 2 - (tc.uR - tc.uL) * (tc.rhoL + tc.rhoR) * (tc.aL + tc.aR) / 8;

    % Two-Shock approximation
    pHat = max(tc.TOL, tc.pPV);
    gL = @(p)sqrt(tc.AL / (p + tc.BL));
    gR = @(p)sqrt(tc.AR / (p + tc.BR));
    tc.pTS = (gL(pHat) * tc.pL + gR(pHat) * tc.pR - (tc.uR - tc.uL)) / (gL(pHat) + gR(pHat));

    % arithmetic mean
    tc.pME = (tc.pL + tc.pR) / 2;
end

% domain of computation
if (pro == 8)
    tc.dm = [0, 1];
elseif (pro == 9)
    tc.dm = [-10, 10];
elseif (pro == 11)
    tc.dm = [-2, 2];
else
    tc.dm = [-1, 1];
end

end



