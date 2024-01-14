% Compute the derivative of velocity with respect to the conserved
% variables
function [du, flag] = getVelocityDeriv(rho, u, v, i, j)

switch i 
    case 1
        switch j
            case 1
                du = -u ./ rho;
                flag = 2;
            case 2
                du = 1 ./ rho;
                flag = 2;
            case 3
                du = zeros(size(rho));
                flag = 0;
            case 4
                du = zeros(size(rho));
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 2
        switch j
            case 1
                du = -v ./ rho;
                flag = 2;
            case 2
                du = zeros(size(rho));
                flag = 0;                
            case 3
                du = 1 ./ rho;
                flag = 2;
            case 4
                du = zeros(size(rho));
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    otherwise
        error('Wrong index')
end

end