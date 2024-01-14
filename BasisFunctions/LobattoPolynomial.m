% Exclude the (1-x.^2)component
% Copyright (C) Xiangyi Meng.

function r = LobattoPolynomial(x, poly_deg, der_order)

switch der_order
    case 0
        switch poly_deg
            case 0
                r = sqrt(3 / 2) * (-2) * ones(size(x));
            case 1
                r = sqrt(5 / 2) * (-2) * x;
            case 2
                r = sqrt(7 / 2) * (5 * x .* x - 1) / (-2);
            case 3
                r = sqrt(9 / 2) * (7 * x .* x - 3) .* x / (-2);
            case 4
                r = sqrt(11 / 2) * ((21 * x .* x - 14) .* x .* x + 1) / (-4);
            case 5
                r = sqrt(13 / 2) * (((33 * x .* x - 30) * x .* x) + 5) .* x / (-4);
            otherwise
                error('Not implemented degree of polynomial')
        end
    case 1
        switch poly_deg
            case 0
                r = zeros(size(x));
            case 1
                r = sqrt(5 / 2) * (-2) * ones(size(x));
            case 2
                r = sqrt(7 / 2) * (-5 * x);
            case 3
                r = sqrt(9 / 2) * (-10.5 * x .* x + 1.5);
            otherwise
                error('Not implemented degree of polynomial')
        end        
    otherwise
        error('Not implemented derivative order of polynomial')
end

end


