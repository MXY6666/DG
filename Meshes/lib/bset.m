% transfers the eight arrays containing 1D boundary data to two 2D arrays (x,y)
function [x, y] = bset(xb_bdry, yb_bdry, xt_bdry, yt_bdry, xl_bdry, yl_bdry, xr_bdry, yr_bdry, m, n)

x = zeros(m + 1, n + 1);
y = zeros(m + 1, n + 1);

x(:, 1) = xb_bdry(:);
y(:, 1) = yb_bdry(:);

x(:, end) = xt_bdry(:);
y(:, end) = yt_bdry(:);

x(1, :) = xl_bdry;
y(1, :) = yl_bdry;

x(end, :) = xr_bdry;
y(end, :) = yr_bdry;

end