% energy flux
function [jE1, jE2] = getEnergyFlux(rho, w1, w2, trz, flag, tc)

tmp = 1 / tc.We .* rho .* trz;
switch flag
    case 1
        jE1 = tmp .* w1;
        jE2 = [];
    case 2
        jE1 = [];
        jE2 = tmp .* w2;
    case 12
        jE1 = tmp .* w1;
        jE2 = tmp .* w2;
    otherwise
        error('Wrong flag in getEnergyFlux')
end

end