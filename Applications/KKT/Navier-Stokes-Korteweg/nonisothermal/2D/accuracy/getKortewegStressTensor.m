% Korteweg stress tensor
function [xi11, xi12, xi22] = getKortewegStressTensor(rho, w1, w2, r, flag, tc)

tmp1 = rho .* r;
tmp2 = 0.5 * (w1.^2 - w2.^2);
xi12 = -1 / tc.We * w1 .* w2;
switch flag
    case 1
        xi11 = (tmp1 - tmp2) / tc.We;       
        xi22 = [];
    case 2
        xi11 = [];
        xi22 = (tmp1 + tmp2) / tc.We;
    case 12
        xi11 = (tmp1 - tmp2) / tc.We;
        xi22 = (tmp1 + tmp2) / tc.We;
    otherwise
        error('Wrong flag in getKortewegStressTensor')
end

end