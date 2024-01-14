% Copyright (C) Xiangyi Meng.

function f = minmod_bar(a, b, c, M, h)
     
isTVB = abs(a) <= M .* h.^2;
f = a .* isTVB;
f = f + minmod(a, b, c) .* ~isTVB;
    
end

