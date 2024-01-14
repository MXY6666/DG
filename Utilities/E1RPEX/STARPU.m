function [P, U] = STARPU(DL, UL, PL, CL, DR, UR, PR, CR, PSCALE,GAMMA)

% *
% C     Purpose: to compute the solution for pressure and
% C              velocity in the Star Region
% *
% C     Declaration of variables
% *
      
      TOLPRE = 1.0E-10;
      NRITER = 20;
% *
% C     Guessed value PSTART is computed
% *

      [PSTART] = GUESSP(DL, UL, PL, CL, DR, UR, PR, CR, GAMMA);

      POLD  = PSTART;
      UDIFF = UR - UL;

%       disp('   Iteration number      Change  ')

      
      iter = 0;
      CHANGE = 1;
      while(CHANGE > TOLPRE && iter < NRITER)


         [FL, FLD] = PREFUN(POLD, DL, PL, CL,GAMMA);
         [FR, FRD] = PREFUN(POLD, DR, PR, CR,GAMMA);

         P      = POLD - (FL + FR + UDIFF)/(FLD + FRD);
         CHANGE = 2.0*abs((P - POLD)/(P + POLD));
%          disp([iter,CHANGE])
         if P<0.0 
             P = TOLPRE;
         end
         POLD  = P;
         iter = iter + 1;
      end
      
      if CHANGE > TOLPRE
          disp('Divergence in Newton-Raphson iteration')
      end


% C     Compute velocity in Star Region

      U = 0.5*(UL + UR + FR - FL);

%       disp('   Pressure        Velocity')
%       disp([P/PSCALE, U])

      return