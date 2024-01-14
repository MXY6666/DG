% Copyright (C) Xiangyi Meng.

function r = LegendrePolynomial(x, poly_deg, der_order)

r = JacobiPolynomial(x, 0, 0, poly_deg, der_order);

% switch der_order
%     case 0
%         switch poly_deg
%             case 0
%                 r = ones(size(x));
%             case 1
%                 r = x;
%             case 2
%                 r = x.^2 - 1 / 3;
%             case 3
%                 r = x.^3 - 3 / 5 * x;
%             case 4
%                 r = x.^4 - 6 / 7 * x.^2 + 3 / 35;
%             case 5
%                 r = x.^5 - 10 / 9 * x.^3 + 5 / 21 * x;
%             case 6
%                 r = x.^6 - 15 / 11 * x.^4 + 5 / 11 * x.^2 - 5 / 231;
%             otherwise
%                 error('Not implemented degree of polynomial')
%         end
%     case 1
%         switch poly_deg
%             case 0
%                 r = zeros(size(x));
%             case 1
%                 r = ones(size(x));
%             case 2
%                 r = 2 * x;
%             case 3
%                 r = 3 * x.^2 - 3 / 5;
%             case 4
%                 r = 4 * x.^3 - 12 / 7 * x;
%             case 5
%                 r = 5 * x.^4 - 10 / 3 * x.^2 + 5 / 21;
%             case 6
%                 r = 6 * x.^5 - 60 / 11 * x.^3 + 10 / 11 * x;
%             otherwise
%                 error('Not implemented degree of polynomial')
%         end        
%     otherwise
%         error('Not implemented derivative order of polynomial')
% end

end
        
        
        
      
        
        