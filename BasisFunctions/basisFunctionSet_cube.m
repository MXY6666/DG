% Copyright (C) Xiangyi Meng.

function r = basisFunctionSet_cube(x, y, z, type, der_order)

deg = mod(type, 10);
if (mod(type, 100) < 10)
    cnt = 1 : deg + 1;
    nb = sum(cnt .* (cnt + 1) / 2);
else
    nb = (deg + 1)^3;
end
np = max(max(numel(x), numel(y)), numel(z));

r = zeros(np, nb);
for j = 1 : nb
    r(:, j) = basisFunction_cube(x(:), y(:), z(:), type, j, der_order);
end

end