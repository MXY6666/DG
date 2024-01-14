% type = 20_ for orthogonal basis functions
% type = 21_ for Lagrange basis functions
% type = 22_ for Bernstein basis functions
% Copyright (C) Xiangyi Meng.

function r = basisFunction_square(x, y, type, index, der_order)

if (size(der_order) ~= 2) 
    error('The argument der_order should be a 2D vector')
end

switch type
    case num2cell(200 : 209)
        deg = mod(type, 10);
        nb = (deg + 1) * (deg + 2) / 2;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end        
        
        [i, j] = getSplitIndex(index);
        r = LegendrePolynomial(x, i, der_order(1)) .* LegendrePolynomial(y, j, der_order(2));
    case num2cell(210 : 214)
        deg = mod(type, 10);
        nb = (deg + 1)^2;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end 
        
        type1 = 110 + deg;
        j = ceil(index / (deg + 1)); 
        if (deg == 1)
            if (mod(index, 3) == 1)
                i = 1;
            else
                i = 2;
            end
        else
            i = index - (j - 1) * (deg + 1);
        end
        r = basisFunction_line(x, type1, i, der_order(1)) .* basisFunction_line(y, type1, j, der_order(2));
    case num2cell(220 : 229)
        deg = mod(type, 10);
        nb = (deg + 1)^2;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end 
        
        type1 = 120 + deg;
        j = ceil(index / (deg + 1)); 
        i = index - (j - 1) * (deg + 1);
        r = basisFunction_line(x, type1, i, der_order(1)) .* basisFunction_line(y, type1, j, der_order(2));        
    otherwise 
        error('Not implemented type of basis functions')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [i, j] = getSplitIndex(index)
    
s = 0;
for n = 1 : 10
    if (s < index && s + n >= index)
        break
    end
    s = s + n;
end
index = index - s;

i = n - index;
j = index - 1;
    
end

% function [i, j] = getSplitIndex(index, deg)
% 
% nb = 0;
% for j = 0 : deg
%     for i = 0 : deg
%        if (i + j <= deg)
%           nb = nb + 1; 
%        end
%        if (nb == index)
%            break
%        end
%     end
%     if (nb == index)
%        break
%     end    
% end
% 
% end







    