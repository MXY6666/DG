% a, b should both be matrices with 2 rows and the same number of columns.
% Copyright (C) Xiangyi Meng.

function c = mycross2(a, b)

c = a(1,:).*b(2,:)-a(2,:).*b(1,:);

end