function  [D,U,P] =  SAMPLE(PM, UM, S,GAMMA,DL, UL, PL, CL, DR, UR, PR, CR)
% *
% C     Purpose: to sample the solution throughout the wave
% C              pattern. Pressure PM and velocity UM in the
% C              Star Region are known. Sampling is performed
% C              in terms of the 'speed' S = X/T. Sampled
% C              values are D, U, P
% *
% C     Input variables : PM, UM, S, /GAMMAS/, /STATES/
% C     Output variables: D, U, P
% *

      G1 = (GAMMA - 1.0)/(2.0*GAMMA);
      G2 = (GAMMA + 1.0)/(2.0*GAMMA);
      G3 = 2.0*GAMMA/(GAMMA - 1.0);
      G4 = 2.0/(GAMMA - 1.0);
      G5 = 2.0/(GAMMA + 1.0);
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0);
      G7 = (GAMMA - 1.0)/2.0;
      G8 = GAMMA - 1.0;
      
      if (S <= UM)
% *
% C        Sampling point lies to the left of the contact
% C        discontinuity
% *
         if (PM <= PL)
% *
% C           Left rarefaction
% *
            SHL = UL - CL;

            if(S <= SHL)
% *
% C              Sampled point is left data state
% *
               D = DL;
               U = UL;
               P = PL;
            else
               CML = CL*(PM/PL)^G1;
               STL = UM - CML;

               if(S > STL)
% *
% C                 Sampled point is Star Left state
% *
                  D = DL*(PM/PL)^(1.0/GAMMA);
                  U = UM;
                  P = PM;
               else
% *
% C                 Sampled point is inside left fan
% *
                  U = G5*(CL + G7*UL + S);
                  C = G5*(CL + G7*(UL - S));
                  D = DL*(C/CL)^G4;
                  P = PL*(C/CL)^G3;
               end
            end
         else
% *
% C           Left shock
% *
            PML = PM/PL;
            SL  = UL - CL*sqrt(G2*PML + G1);

            if(S <= SL)
% *
% C              Sampled point is left data state
% *
               D = DL;
               U = UL;
               P = PL;

            else
% *
% C              Sampled point is Star Left state
% *
               D = DL*(PML + G6)/(PML*G6 + 1.0);
               U = UM;
               P = PM;
            end
         end
      else
% *
% C        Sampling point lies to the right of the contact discontinuity
% *
         if(PM>PR)
% *
% C           Right shock
% *
            PMR = PM/PR;
            SR  = UR + CR*sqrt(G2*PMR + G1);

            if(S>=SR)
% *
% C              Sampled point is right data state
% *
               D = DR;
               U = UR;
               P = PR;
            else
% *
% C              Sampled point is Star Right state
% *
               D = DR*(PMR + G6)/(PMR*G6 + 1.0);
               U = UM;
               P = PM;
            end
         else
% *
% C           Right rarefaction
% *
            SHR = UR + CR;

            if(S>=SHR)
% *
% C              Sampled point is right data state
% *
               D = DR;
               U = UR;
               P = PR;
            else
               CMR = CR*(PM/PR)^G1;
               STR = UM + CMR;

               if(S<=STR)
% *
% C                 Sampled point is Star Right state
% *
                  D = DR*(PM/PR)^(1.0/GAMMA);
                  U = UM;
                  P = PM;
               else
% *
% C                 Sampled point is inside left fan
% *
                  U = G5*(-CR + G7*UR + S);
                  C = G5*(CR - G7*(UR - S));
                  D = DR*(C/CR)^G4;
                  P = PR*(C/CR)^G3;
               end
            end
         end
      end

      return