% Compute the derivative of velocity with respect to the conserved
% variables
function [du, flag] = getVelocityDeriv(rho, u, i)

switch i
    case 1
        du = -u ./ rho;
        flag = 2;
    case 2
        du = 1 ./ rho;
        flag = 2;
    case 3
        du = zeros(size(rho));
        flag = 0;
    otherwise
        error('Wrong index')
end

end