% Copyright (C) Xiangyi Meng.

function r = basisFunctionSet_square(x, y, type, der_order)

deg = mod(type, 10);
if (mod(type, 100) < 10)
    nb = (deg + 1) * (deg + 2) / 2;
else
    nb = (deg + 1)^2;
end
np = max(numel(x), numel(y));

r = zeros(np, nb);
for j = 1 : nb
    r(:, j) = basisFunction_square(x(:), y(:), type, j, der_order);
end

end