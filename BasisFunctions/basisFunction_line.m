% type = 10_ for orthogonal basis functions
% type = 11_ for Lagrange basis functions
% type = 12_ for Bernstein basis functions
% Copyright (C) Xiangyi Meng.

function r = basisFunction_line(x, type, index, der_order)

switch type
    case num2cell(100 : 109)
        deg = mod(type, 10);
        nb = deg + 1;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end
        
        r = LegendrePolynomial(x, index - 1, der_order);
    case 110
        if (index ~= 1)
            error('There is only one basis functions for this type')
        end
        
        switch der_order
            case 0
                r = ones(size(x));
            case 1
                r = zeros(size(x));
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    case 111
        switch der_order
            case 0
                switch index
                    case 1
                        r = (1 - x) / 2;
                    case 2
                        r = (1 + x) / 2;
                end       
            case 1
                switch index
                    case 1
                        r = -1 / 2 * ones(size(x));
                    case 2
                        r = 1 / 2 * ones(size(x));
                end                
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    case 112
        switch der_order
            case 0
                switch index
                    case 1
                        r = x .* (x - 1) / 2;
                    case 2
                        r = (1 + x) .* (1 - x);
                    case 3
                        r = x .* (x + 1) / 2;
                end    
            case 1
                switch index
                    case 1
                        r = x - 1 / 2;
                    case 2
                        r = -2 * x;
                    case 3
                        r = x + 1 / 2;
                end                
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    case 113
        switch der_order
            case 0
                switch index
                    case 1
                        r = -9 / 16 * (x + 1 / 3) .* (x - 1 / 3) .* (x - 1);
                    case 2
                        r = 27 / 16 * (x + 1) .* (x - 1 / 3) .* (x - 1);
                    case 3
                        r = -27 / 16 * (x + 1) .* (x + 1 / 3) .* (x - 1);
                    case 4
                        r = 9 / 16 * (x + 1) .* (x + 1 / 3) .* (x - 1 / 3);
                end  
            case 1
                switch index
                    case 1
                        r = -27 / 16 * x.^2 + 9 / 8 * x  + 1 / 16;
                    case 2
                        r = 81 / 16 * x.^2 - 9 / 8 * x - 27 / 16;
                    case 3
                        r = -81 / 16 * x.^2 - 9 / 8 * x + 27 / 16;
                    case 4
                        r = 27 / 16 * x.^2 + 9 / 8 * x - 1 / 16;
                end                
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    case 114
        switch der_order
            case 0
                switch index
                    case 1
                        r = 2 / 3 * (x + 1 / 2) .* x .* (x - 1 / 2) .* (x - 1);
                    case 2
                        r = -8 / 3 * (x + 1) .* x .* (x - 1 / 2) .* (x - 1);
                    case 3
                        r = 4 * (x + 1) .* (x + 1 / 2) .* (x - 1 / 2) .* (x - 1);
                    case 4
                        r = -8 / 3 * (x + 1) .* (x + 1 / 2) .* x .* (x - 1);
                    case 5
                        r = 2 / 3 * (x + 1) .* (x + 1 / 2) .* x .* (x - 1 / 2);
                end
            case 1
                switch index
                    case 1
                        r = (8 * x.^3) / 3 - 2 * x.^2 - x / 3 + 1 / 6;
                    case 2
                        r = -(32 * x.^3) / 3 + 4 * x.^2 + (16 * x) / 3 - 4 / 3;
                    case 3
                        r = 2 * x .* (8 * x.^2 - 5);
                    case 4
                        r = -(32 * x.^3) / 3 - 4 * x.^2 + (16 * x) / 3 + 4 / 3;
                    case 5
                        r = (8 * x.^3) / 3 + 2 * x.^2 - x / 3 - 1 / 6;
                end                
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    case num2cell(120 : 129)
        deg = mod(type, 10);
        nb = deg + 1;
        if all(index ~= 1 : nb)
            error(['There is only ', num2str(nb), ' basis functions for this type'])
        end
        
        switch der_order
            case 0
                r = BernsteinPolynomial(0.5 * (x + 1), deg, index - 1, 0);
            case 1
                r = 0.5 * BernsteinPolynomial(0.5 * (x + 1), deg, index - 1, 1);
            otherwise
                error('Not implemented derivative order of basis funtions')
        end
    otherwise
        error('Not implemented type of basis functions')
end

end