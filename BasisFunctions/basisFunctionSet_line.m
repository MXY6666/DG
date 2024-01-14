% Copyright (C) Xiangyi Meng.

function r = basisFunctionSet_line(x, type, der_order)

deg = mod(type, 10);
nb = deg + 1;
np = numel(x);

r = zeros(np, nb);
for j = 1 : nb
    r(:, j) = basisFunction_line(x(:), type, j, der_order);
end

end