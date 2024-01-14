% Compute the derivative of -kappa * T with respect to the conserved
% variables
function [dsT, flag] = getScaledTemperatureDeriv(rho, u, v, E, i, tc)

switch i
    case 1
        dsT = tc.gamma * tc.mu / tc.Pr * (E ./ rho - u.^2 - v.^2) ./ rho;
    case 2
        dsT = tc.gamma * tc.mu / tc.Pr * u ./ rho;
    case 3
        dsT = tc.gamma * tc.mu / tc.Pr * v ./ rho;
    case 4
        dsT = -tc.gamma * tc.mu / tc.Pr ./ rho;
    otherwise
        error('Wrong index')
end
flag = 2;

end