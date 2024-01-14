% Copyright (C) Xiangyi Meng.

function f = minmod(a, b, c)

f = min(abs(a), abs(b));
f = min(f, abs(c));

f = sign(a) .* f .* (sign(a) == sign(b) & sign(a) == sign(c));

end

