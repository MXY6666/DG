function [PM] = GUESSP(DL, UL, PL, CL, DR, UR, PR, CR, GAMMA)
% *
% C     Purpose: to provide a guessed value for pressure
% C              PM in the Star Region. The choice is made
% C              according to adaptive Riemann solver using
% C              the PVRS, TRRS and TSRS approximate
% C              Riemann solvers. See Sect. 9.5 of Chapt. 9
% C              of Ref. 1
% *


      
      G1 = (GAMMA - 1.0)/(2.0*GAMMA);
      G2 = (GAMMA + 1.0)/(2.0*GAMMA);
      G3 = 2.0*GAMMA/(GAMMA - 1.0);
      G4 = 2.0/(GAMMA - 1.0);
      G5 = 2.0/(GAMMA + 1.0);
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0);
      G7 = (GAMMA - 1.0)/2.0;
      G8 = GAMMA - 1.0;

      QUSER = 2.0;

% C     Compute guess pressure from PVRS Riemann solver

      CUP  = 0.25*(DL + DR)*(CL + CR);
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP;
      PPV  = max(0.0, PPV);
      PMIN = min(PL,  PR);
      PMAX = max(PL,  PR);
      QMAX = PMAX/PMIN;

      if ((QMAX <= QUSER) &&  (PMIN<=PPV && PPV<=PMAX))

% C        Select PVRS Riemann solver
         PM = PPV;
         
      else
         if (PPV<PMIN)

% C           Select Two-Rarefaction Riemann solver

            PQ  = (PL/PR)^G1;
            UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
            PTL = 1.0 + G7*(UL - UM)/CL;
            PTR = 1.0 + G7*(UM - UR)/CR;
            PM  = 0.5*(PL*PTL^G3 + PR*PTR^G3);
         else

% C           Select Two-Shock Riemann solver with PVRS as estimate

            GEL = sqrt((G5/DL)/(G6*PL + PPV));
            GER = sqrt((G5/DR)/(G6*PR + PPV));
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
         end
      end

      return