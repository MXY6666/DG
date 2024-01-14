% performs 2D transfinite interpolation on the eight boundary
% arrays to obtain the grid (x, y)
function [x, y] = tfi(x, y, xb, yb, xt, yt, xl, yl, xr, yr, m, n, dx, dy)

for j = 2 : n
    s = (j - 1) * dy;
    for i = 2 : m
        r = (i - 1) * dx;
        x(i, j) = (1 - s) * xb(i) + s * xt(i) + (1 - r) * xl(j) + r * xr(j) ...
                - (r * s * xr(n + 1) + (1 - r) * s * xl(n + 1) ...
                + r * (1 - s) * xr(1) + (1 - r) * (1 - s) * xl(1));
        y(i, j) = (1 - s) * yb(i) + s * yt(i) + (1 - r) * yl(j) + r * yr(j) ...
                - (r * s * yr(n + 1) + (1 - r) * s * yl(n + 1) ...
                + r * (1 - s) * yr(1) + (1 - r) * (1 - s) * yl(1));
    end
end

end
    