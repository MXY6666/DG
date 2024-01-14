function [x, y] = tensor(m, n, idmn, a, it, imth)

if (nargin < 2)
    m = input('Input number of elements in xi direction : ');
    n = input('Input number of elements in eta direction : ');
end

% grid size
[dx, dy, tdx, tdy, dxx, dxy, dyy] = discrtz(m, n);

% Select physical region
if (nargin < 4)  || isempty(idmn) || isempty(a)
    disp('Physical Domain: ')
    [idmn, a1, b1, a2, b2, a3, b3, a4, b4] = select;
else
    a1 = a(1, 1); b1 = a(2, 1);
    a2 = a(1, 2); b2 = a(2, 2);
    a3 = a(1, 3); b3 = a(2, 3);
    a4 = a(1, 4); b4 = a(2, 4);
end

% discretizations for four boundaries of reference domain
if (nargin < 5)
    it = [];
end
[tb, tt, tl, tr] = ts(m, n, it);

% discretizations for 'four' boundaries of physical domain
[xb_bdry, yb_bdry, xt_bdry, yt_bdry, xl_bdry, yl_bdry, xr_bdry, yr_bdry]...
= b_parm(m, n, idmn, a1, b1, a2, b2, a3, b3, a4, b4, tb, tt, tl, tr);

% Insert the eight boundary vectors to the two discretization matrices of
% the whole domain
[x, y] = bset(xb_bdry, yb_bdry, xt_bdry, yt_bdry, xl_bdry, yl_bdry, xr_bdry, yr_bdry, m, n);

% Generate initial grid
[x, y] = tfi(x, y, xb_bdry, yb_bdry, xt_bdry, yt_bdry, xl_bdry, yl_bdry, xr_bdry,yr_bdry, m, n, dx, dy);

if (nargin < 6)  || isempty(imth)
    disp(' 1) AO ')
    disp(' 2) TTM ')
    disp(' 3) Length ')
    disp(' 4) Area ')
    disp(' 5) 0.9 Area + 0.1 Smo. ')
    disp(' 6) g12^2 ')
    disp(' 7) Liao ')
    disp(' 8) Mod-Liao ')
    disp(' 9) Scaled-Laplace ')
    disp('10) g12*g12/g11/g22 ')
    disp('11) g11+g22 -rg() ')
    disp('12) AO**2 ')
    disp('13) Guest Functional ')
    imth = input('Input Method: ');
end

% Solve pde
[rhsx, rhsy] = sset(m, n, idmn, tb, tt, tl, tr);

icntmax = 10000;
ncntmax = 6;
tol = 1.e-05;
omg = 1;

icnt = 0;
xerr = 1;
while (xerr > tol)
    icnt = icnt + 1;
    [x, y, xerr] = sor(x, y, rhsx, rhsy, dxx, dxy, dyy, tdx, tdy, m, n, ncntmax, omg, imth);
    fprintf('icnt = %d,   xerr = %f,   tol = %f\n', icnt, xerr, tol)
    
    if (icnt > icntmax)
        disp('Outer Iteration Count Exceeded')
        break;
    end
end
if (icnt ~= icntmax)
    fprintf('tol satisfied after icnt = %d\n', icnt)
end

% final grid
fprintf('\n\n                       Final Grid\n')
grd_met(x, y, m, n)

end


