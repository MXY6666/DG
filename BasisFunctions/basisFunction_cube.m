% type = 30_ for orthogonal basis functions
% type = 31_ for Lagrange basis functions
% type = 32_ for Bernstein basis functions
% Copyright (C) Xiangyi Meng.

function r = basisFunction_cube(x, y, z, type, index, der_order)

if (size(der_order) ~= 3) 
    error('The argument der_order should be a 3D vector')
end

switch type
    case num2cell(300 : 309)
        deg = mod(type, 10);
        cnt = 1 : deg + 1;
        nb = sum(cnt .* (cnt + 1) / 2);
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end   
        
        [i, j, k] = getSplitIndex(index);
        r = LegendrePolynomial(x, i, der_order(1)) .* LegendrePolynomial(y, j, der_order(2)) .* LegendrePolynomial(z, k, der_order(3));
    case num2cell(310 : 314)
        deg = mod(type, 10);
        nb = (deg + 1)^3;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end 
        
        type1 = 110 + deg;
        k = ceil(index / (deg + 1)^2);
        j = ceil((index - (k - 1) * (deg + 1)^2) / (deg + 1)); 
        if (deg == 1)
            if (mod(index - (k - 1) * (deg + 1)^2, 3) == 1)
                i = 1;
            else
                i = 2;
            end
        else
            i = index - (k - 1) * (deg + 1)^2 - (j - 1) * (deg + 1);
        end                
        r = basisFunction_line(x, type1, i, der_order(1)) .* basisFunction_line(y, type1, j, der_order(2)) .* basisFunction_line(z, type1, k, der_order(3));
    case num2cell(320 : 329)
        deg = mod(type, 10);
        nb = (deg + 1)^3;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end 
        
        type1 = 120 + deg;
        k = ceil(index / (deg + 1)^2);
        j = ceil((index - (k - 1) * (deg + 1)^2) / (deg + 1)); 
        i = index - (k - 1) * (deg + 1)^2 - (j - 1) * (deg + 1);            
        r = basisFunction_line(x, type1, i, der_order(1)) .* basisFunction_line(y, type1, j, der_order(2)) .* basisFunction_line(z, type1, k, der_order(3));
    otherwise 
        error('Not implemented type of basis functions')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [i, j, k] = getSplitIndex(index)
    
s = 0;
for n = 1 : 10
    tmp = n * (n + 1) / 2;
    if (s < index && s + tmp >= index)
        break
    end
    s = s + tmp;
end
index = index - s;

s = 0;
for l = 1 : n
    if (s < index && s + l >= index)
        break
    end
    s = s + l;
end
index = index - s;

i = n - l;
j = l - index;
k = index - 1;

end

end



