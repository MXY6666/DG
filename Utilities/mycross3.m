% Simplify the internal function 'cross'
% a, b should both be matrices with 3 rows and the same number of columns.
% Copyright (C) Xiangyi Meng.

function c = mycross3(a, b)

c(3,:) = a(1,:).*b(2,:)-a(2,:).*b(1,:);
c(1,:) = a(2,:).*b(3,:)-a(3,:).*b(2,:);
c(2,:) = a(3,:).*b(1,:)-a(1,:).*b(3,:);

end