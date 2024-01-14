function [CFL, CFLMax] = setCFLNumber(tc)

switch tc.pro 
    % Sod problem, left rarefaction, middle contact discontinuity, right shock
    case 1
        CFL    = 0.1; 
        CFLMax = 1;
    % Lax problem, left rarefaction, middle contact discontinuity, right shock
    case 2
        CFL    = 0.1; 
        CFLMax = 1;
                   
    % left rarefaction, middle contact discontinuity, right shock
    case 3
        CFL    = 0.1; 
        CFLMax = 1;  
        
    % left shock, middle contact discontinuity, right rarefaction
    case 4
        CFL    = 0.1; 
        CFLMax = 1;   
       
    % left shock, middle contact discontinuity, right shock
    case 5
        CFL    = 0.1; 
        CFLMax = 1; 
         
    % 123 problem, two strong rarefactions and a stationary contact discontinuity
    case 6
        CFL    = 0.1; 
        CFLMax = 1;
        
    % 123 problem, two strong rarefactions and a stationary contact discontinuity
    case 7
        CFL    = 0.1; 
        CFLMax = 1; 
        
    % Leblanc problem
    case 8
        CFL    = 1.e-4; 
        CFLMax = 1;
        
    % Leblanc problem
    case 9
        CFL    = 1.e-10; 
        CFLMax = 1;
        
    % Sedov blast wave problem
    case 10
        CFL    = 1.e-5; 
        CFLMax = 1;
        
    % Sedov blast wave problem
    case 11 
        CFL    = 1.e-10; 
        CFLMax = 1;
    otherwise
        error('Not implemented problem')
end 

end



