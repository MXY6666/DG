function  [F,FD] = PREFUN(P,DK,PK,CK,GAMMA)
% *
% C     Purpose: to evaluate the pressure functions
% C              FL and FR in exact Riemann solver
% C              and their first derivatives
% *
     
      G1 = (GAMMA - 1.0)/(2.0*GAMMA);
      G2 = (GAMMA + 1.0)/(2.0*GAMMA);
      G3 = 2.0*GAMMA/(GAMMA - 1.0);
      G4 = 2.0/(GAMMA - 1.0);
      G5 = 2.0/(GAMMA + 1.0);
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0);
      G7 = (GAMMA - 1.0)/2.0;
      G8 = GAMMA - 1.0;

      if (P <= PK)
% *
% C        Rarefaction wave
% *
         PRATIO = P/PK;
         F    = G4*CK*(PRATIO^G1 - 1.0);
         FD   = (1.0/(DK*CK))*PRATIO^(-G2);
      else
% *
% C        Shock wave
% *
         AK  = G5/DK;
         BK  = G6*PK;
         QRT = sqrt(AK/(BK + P));
         F   = (P - PK)*QRT;
         FD  = (1.0 - 0.5*(P - PK)/(BK + P))*QRT;
      end

      return