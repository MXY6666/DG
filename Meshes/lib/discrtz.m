%   m cells in xi-direction, n cells in eta-direction
function [dx, dy, tdx, tdy, dxx, dxy, dyy] = discrtz(m, n)
 
dx = 1 / m;
dy = 1 / n;
tdx = 2 * dx;
tdy = 2 * dy;
dxx = dx * dx;
dxy = 4 * dx * dy;
dyy = dy  * dy;
   
end