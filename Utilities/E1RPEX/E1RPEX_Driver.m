% Exact Riemann solver driver

% %C---------------------------------------------------------------------C

% DOMLEN = 1.0;     %  DOMLEN : Domain length,               TEST 1 (Modified Sod)
% DIAPH1 = 0.3;     %  DIAPH1 : Position of diaphragm   
% CELLS = 1000;    %  CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;     %  GAMMA  : Ratio of specific heats
% TIMEOU = 0.20;    %  TIMEOU : Output time
% DL = 1.0;     %  DL     : Initial density  on left  section of tube
% UL = 0.75;    %  UL     : Initial velocity on left  section of tube
% PL = 1.0;     %  PL     : Initial pressure on left  section of tube
% DR = 0.125;   %  DR     : Initial density  on right section of tube
% UR = 0.0;     %  UR     : Initial velocity on right section of tube
% PR = 0.1;     %  PR     : Initial pressure on right section of tube
% PSCALE = 1.0;     %  PSCALE : Normalising factor for pressure and energy

% %C---------------------------------------------------------------------C
% TEST 2 (123 problem)
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.5;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.15;        % TIMEOU : Output time
% DL = 1.0;         % DL     : Initial density  on left  section of tube
% UL = -2.0;        % UL     : Initial velocity on left  section of tube
% PL = 0.4;         % PL     : Initial pressure on left  section of tube
% DR = 1.0;         % DR     : Initial density  on right section of tube
% UR = 2.0;         % UR     : Initial velocity on right section of tube
% PR = 0.4;         % PR     : Initial pressure on right section of tube
% PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy
% 
% %C---------------------------------------------------------------------C
% TEST 3 (Left Woodward & Colella)
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.5;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.012;        % TIMEOU : Output time
% DL = 1.0;         % DL     : Initial density  on left  section of tube
% UL = 0.0;        % UL     : Initial velocity on left  section of tube
% PL = 1000.0;         % PL     : Initial pressure on left  section of tube
% DR = 1.0;         % DR     : Initial density  on right section of tube
% UR = 0.0;         % UR     : Initial velocity on right section of tube
% PR = 0.01;         % PR     : Initial pressure on right section of tube
% PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy
% 
% %C---------------------------------------------------------------------C
% TEST 4 (Collision of 2 shocks)
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.4;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.035;        % TIMEOU : Output time
% DL = 5.99924;         % DL     : Initial density  on left  section of tube
% UL = 19.5975;        % UL     : Initial velocity on left  section of tube
% PL = 460.894;         % PL     : Initial pressure on left  section of tube
% DR = 5.99242;         % DR     : Initial density  on right section of tube
% UR = -6.19633;         % UR     : Initial velocity on right section of tube
% PR = 46.0950;         % PR     : Initial pressure on right section of tube
% PSCALE = 1.0;     % PSCALE : Normalising factor for pressure and energy
% 
% %C---------------------------------------------------------------------C
% TEST 5 (Stationary contact)
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.8;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.012;        % TIMEOU : Output time
% DL = 1.0;         % DL     : Initial density  on left  section of tube
% UL = -19.59745;   % UL     : Initial velocity on left  section of tube
% PL = 1000.0;      % PL     : Initial pressure on left  section of tube
% DR = 1.0;         % DR     : Initial density  on right section of tube
% UR = -19.59745   % UR     : Initial velocity on right section of tube
% PR = 0.01        % PR     : Initial pressure on right section of tube
% PSCALE = 1.0E+00     % PSCALE : Normalising factor for pressure and energy

% %C---------------------------------------------------------------------C
% TEST 5 (Normal moving shock) (Steady) ?
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.5;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.3;        % TIMEOU : Output time
% 
% 
% 
% 
% DL = 1.4;         % DL     : Initial density  on left  section of tube
% UL = 10;   % UL     : Initial velocity on left  section of tube
% PL = 1;      % PL     : Initial pressure on left  section of tube
% 
% AL = sqrt(GAMMA*PL/DL);  % sound speed on the left
% ML = UL/AL;
% 
% 
% Cd = (GAMMA+1)*ML^2/( (GAMMA-1)*ML^2+2 );
% Cp = (2*GAMMA*ML^2 - (GAMMA-1) )/(GAMMA+1);
% 
% DR = DL*Cd;         % DR     : Initial density  on right section of tube
% PR = PL*Cp;        % PR     : Initial pressure on right section of tube
% AR = sqrt(GAMMA*PR/DR);  % sound speed on the right
% MR2 = ( (GAMMA-1)*ML^2+2 )/(2*GAMMA*ML^2-(GAMMA-1));
% MR = sqrt(MR2);  %Mach number on the right
% UR = MR*AR;      % UR     : Initial velocity on right section of tube
% 
% PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy

% % %C---------------------------------------------------------------------C
% % TEST 6 (moving shock) (Steady) ?
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.5;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.02;        % TIMEOU : Output time
% 
% W = 5.09; % relative shock speed
% 
% DR = 1.4;
% UR = 0;
% PR = 1;
% AR = sqrt(GAMMA*PR/DR);
% MR = UR/AR;
% 
% Mx = W/AR;
% AL = AR * sqrt( 1+ 2*(GAMMA-1)/(GAMMA+1)^2*( GAMMA*Mx^2-1/Mx^2-(GAMMA-1) ) );
% 
% DL = DR/(1-2/(GAMMA+1)*(1-1/Mx^2));
% PL = PR* (1+2*GAMMA/(GAMMA+1)*(Mx^2-1));
% UL = Mx*(1-((GAMMA-1)*Mx^2 +2)/((GAMMA+1)*Mx^2) )*AR;  
% 
% 
% 
% PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy
% % 
% %C---------------------------------------------------------------------C
% TEST 5 (Normal moving shock) inverse (Steady) not correct ?
% DOMLEN = 1.0;         % DOMLEN : Domain length,
% DIAPH1 =0.5;         % DIAPH1 : Position of diaphragm 
% CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
% GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
% TIMEOU = 0.3;        % TIMEOU : Output time
% 

% 
% 
% 
% DR = 1.4;         % DL     : Initial density  on left  section of tube
% UR = 0;   % UL     : Initial velocity on left  section of tube
% PR = 1;      % PL     : Initial pressure on left  section of tube
% 
% AR = sqrt(GAMMA*PL/DL);  % sound speed on the left
% MR = 0;
% 
% ML = 10;
% 
% 
% Cd = (GAMMA+1)*ML^2/( (GAMMA-1)*ML^2+2 );
% Cp = (2*GAMMA*ML^2 - (GAMMA-1) )/(GAMMA+1);
% 
% DL = DR/Cd;         % DR     : Initial density  on right section of tube
% PL = PR/Cp;        % PR     : Initial pressure on right section of tube
% AL = sqrt(GAMMA*PL/DL);  % sound speed on the right
% %MR2 = ( (GAMMA-1)*ML^2+2 )/(2*GAMMA*ML^2-(GAMMA-1));
% %MR = sqrt(MR2);  %Mach number on the right
% UL = ML*AL;      % UR     : Initial velocity on right section of tube
% 
% PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy

% %C---------------------------------------------------------------------C
% TEST (LeBlanc)
DOMLEN = 1.0;         % DOMLEN : Domain length,
DIAPH1 =0.3;         % DIAPH1 : Position of diaphragm 
CELLS = 1000;        % CELLS  : Number of cells in evaluating exact solution
GAMMA = 1.4;         % GAMMA  : Ratio of specific heats
TIMEOU = 0.12;        % TIMEOU : Output time
DL = 1.0E+04;         % DL     : Initial density  on left  section of tube
UL = 0;   % UL     : Initial velocity on left  section of tube
PL = 1.0E+04;      % PL     : Initial pressure on left  section of tube
DR = 1.0;         % DR     : Initial density  on right section of tube
UR = 0;   % UR     : Initial velocity on right section of tube
PR = 1.0;        % PR     : Initial pressure on right section of tube
PSCALE = 1.0E+00;     % PSCALE : Normalising factor for pressure and energy




[XS,DS,US,PS,ES] = E1RPEX(DOMLEN, DIAPH1, CELLS, GAMMA, TIMEOU, DL, UL, PL, DR, UR, PR, PSCALE);

h=1/1000;
xmesh=h/2:h:1-h/2;
hold on;
semilogy(xmesh,DS,'r')
%hold on;semilogy(xmesh,abs(US)+sqrt(1.4*PS./DS),'r')
