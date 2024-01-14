% Given the physical elements, mapping the physical points in them to the 
% points in the reference element.
% Vertex coordinates of each element are stored in each column of xx, yy,
% respectively. Physical point coordinates in each element are stored in
% each column of x, y, respectively.
% Copyright (C) Xiangyi Meng.

function [xi, eta] = phys2Ref_quad(x, y, xx, yy)

a0 = 0.25 * sum(xx);
a1 = 0.25 * (-xx(1, :) + xx(2, :) + xx(3, :) - xx(4, :));
a2 = 0.25 * (-xx(1, :) - xx(2, :) + xx(3, :) + xx(4, :));
a3 = 0.25 * (+xx(1, :) - xx(2, :) + xx(3, :) - xx(4, :));

b0 = 0.25 * sum(yy);
b1 = 0.25 * (-yy(1, :) + yy(2, :) + yy(3, :) - yy(4, :));
b2 = 0.25 * (-yy(1, :) - yy(2, :) + yy(3, :) + yy(4, :));
b3 = 0.25 * (+yy(1, :) - yy(2, :) + yy(3, :) - yy(4, :));

xIsDeg = (a2 == 0) & (a3 == 0);
yIsDeg = (b1 == 0) & (b3 == 0);
ind1 = find(xIsDeg & yIsDeg);
ind2 = find(xIsDeg & ~yIsDeg);
ind3 = find(~xIsDeg & yIsDeg);
ind4 = find(~xIsDeg & ~yIsDeg);

xi = zeros(size(x));
eta = zeros(size(y));
if ~isempty(ind1)
    xi(:, ind1)  = (x(:, ind1) - a0(ind1)) ./ a1(ind1);
    eta(:, ind1) = (y(:, ind1) - b0(ind1)) ./ b2(ind1);
end
if ~isempty(ind2)
    temp = x(:, ind2) - a0(ind2);
    xi(:, ind2)  = temp ./ a1(ind2);
    eta(:, ind2) = (a1(ind2) .* (y(:, ind2) - b0(ind2)) - b1(ind2) .* temp) ./ (a1(ind2) .* b2(ind2) + b3(ind2) .* temp);
end
if ~isempty(ind3)
    temp = y(:, ind3) - b0(ind3);
    xi(:, ind3)  = -(a2(ind3) .* temp - b2(ind3) .* (x(:, ind3) - a0(ind3))) ./ (a1(ind3) .* b2(ind3) + a3(ind3) .* temp);
    eta(:, ind3) = temp ./ b2(ind3);
end
if ~isempty(ind4)
    temp1 = a2(ind4) .* b1(ind4) - a1(ind4) .* b2(ind4);
    temp2 = a3(ind4) .* b0(ind4) - a0(ind4) .* b3(ind4) + b3(ind4) .* x(:, ind4) - a3(ind4) .* y(:, ind4);
    C1 = 2 * (a3(ind4) .* b1(ind4) - a1(ind4) .* b3(ind4));
    C2 = 2 * (a3(ind4) .* b2(ind4) - a2(ind4) .* b3(ind4));
    A1 = temp1 + temp2;
    A2 = temp1 - temp2;
    B = sqrt(A1.^2 + 2 * C1 .* (b2(ind4) .* (a0(ind4) - x(:, ind4)) + a2(ind4) .* (y(:, ind4) - b0(ind4))));

    xi(:, ind4)  = -(A1 + B) ./ C1;
    eta(:, ind4) = (A2 + B) ./ C2;    
end

end



