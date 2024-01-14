function src = getSourceTerm(x, y, t, tc)

rho    = tc.rho(x, y, t);
rhox   = tc.rho_x(x, y, t);
rhoy   = tc.rho_y(x, y, t);
rhot   = tc.rho_t(x, y, t);
rhoxx  = tc.rho_xx(x, y, t);
rhoyy  = tc.rho_yy(x, y, t);
rhoxy  = tc.rho_xy(x, y, t);
rhoxt  = tc.rho_xt(x, y, t);
rhoyt  = tc.rho_yt(x, y, t);
rhoxxx = tc.rho_xxx(x, y, t);
rhoyyy = tc.rho_yyy(x, y, t);
rhoxxy = tc.rho_xxy(x, y, t);
rhoxyy = tc.rho_xyy(x, y, t);

u   = tc.u(x, y, t);
ux  = tc.u_x(x, y, t);
uy  = tc.u_y(x, y, t);
ut  = tc.u_t(x, y, t);
uxx = tc.u_xx(x, y, t);
uxy = tc.u_xy(x, y, t);
uyy = tc.u_yy(x, y, t);

v   = tc.v(x, y, t);
vx  = tc.v_x(x, y, t);
vy  = tc.v_y(x, y, t);
vt  = tc.v_t(x, y, t);
vxx = tc.v_xx(x, y, t);
vxy = tc.v_xy(x, y, t);
vyy = tc.v_yy(x, y, t);

T   = tc.T(x, y, t);
Tx  = tc.T_x(x, y, t);
Ty  = tc.T_y(x, y, t);
Tt  = tc.T_t(x, y, t);
Txx = tc.T_xx(x, y, t);
Tyy = tc.T_yy(x, y, t);

m1     = rho .* u;
m2     = rho .* v;
u2     = u.^2;
uv     = u .* v;
v2     = v.^2;
q2     = u2 + v2;
r      = rhoxx + rhoyy;
d      = ux + vy;
tmp1   = 8 / 3 * tc.Cv * T - 3 * rho + 0.5 * q2;
E      = rho  .* tmp1 + 0.5 / tc.We * (rhox.^2 + rhoy.^2);
Ex     = rhox .* tmp1 + rho .* (8 / 3 * tc.Cv * Tx - 3 * rhox + u .* ux + v .* vx) + (rhox .* rhoxx + rhoy .* rhoxy) / tc.We;
Ey     = rhoy .* tmp1 + rho .* (8 / 3 * tc.Cv * Ty - 3 * rhoy + u .* uy + v .* vy) + (rhox .* rhoxy + rhoy .* rhoyy) / tc.We;
Et     = rhot .* tmp1 + rho .* (8 / 3 * tc.Cv * Tt - 3 * rhot + u .* ut + v .* vt) + (rhox .* rhoxt + rhoy .* rhoyt) / tc.We;
tmp1   = 8 * rho ./ (3 - rho);
p      = tmp1 .* T - 3 * rho.^2;
px     = 24 * rhox ./ (3 - rho).^2 .* T + tmp1 .* Tx - 6 * rho .* rhox;
py     = 24 * rhoy ./ (3 - rho).^2 .* T + tmp1 .* Ty - 6 * rho .* rhoy;
tau11  = (4 / 3 * ux - 2 / 3 * vy) / tc.Re;
tau22  = (4 / 3 * vy - 2 / 3 * ux) / tc.Re;
tau12  = (uy + vx) / tc.Re;
tau11x = (4 / 3 * uxx - 2 / 3 * vxy) / tc.Re;
tau22y = (4 / 3 * vyy - 2 / 3 * uxy) / tc.Re;
tau12x = (uxy + vxx) / tc.Re;
tau12y = (uyy + vxy) / tc.Re;
tmp1   = rho .* r;
tmp2   = 0.5 * (rhox.^2 - rhoy.^2);
xi11   = (tmp1 - tmp2) / tc.We;
xi22   = (tmp1 + tmp2) / tc.We;
xi12   = -1 / tc.We * rhox .* rhoy;
xi11x  = (rhox .* r + rho .* (rhoxxx + rhoxyy) - (rhox .* rhoxx - rhoy .* rhoxy)) / tc.We;
xi22y  = (rhoy .* r + rho .* (rhoxxy + rhoyyy) - (rhox .* rhoxy - rhoy .* rhoyy)) / tc.We;
xi12x  = -1 / tc.We * (rhoxx .* rhoy + rhox .* rhoxy);
xi12y  = -1 / tc.We * (rhoxy .* rhoy + rhox .* rhoyy);
q1x    = -(8 * tc.Cv) / (3 * tc.We * tc.Pr) * Txx;
q2y    = -(8 * tc.Cv) / (3 * tc.We * tc.Pr) * Tyy;
jE1x   = 1 / tc.We * (rhox.^2 .* d + rho .* rhoxx .* d + rho .* rhox .* (uxx + vxy));
jE2y   = 1 / tc.We * (rhoy.^2 .* d + rho .* rhoyy .* d + rho .* rhoy .* (uxy + vyy));

src = [rhot + rho .* d + rhox .* u + rhoy .* v; ...
       rhot .* u + rho .* ut + rhox .* u2 + rhoy .* uv + m1 .* (ux + d) + m2 .* uy + px - tau11x - tau12y - xi11x - xi12y; ...
       rhot .* v + rho .* vt + rhox .* uv + rhoy .* v2 + m2 .* (vy + d) + m1 .* vx + py - tau12x - tau22y - xi12x - xi22y; ...
       Et + u .* (Ex + px - tau11x - tau12y - xi11x - xi12y) + v .* (Ey + py - tau12x - tau22y - xi12x - xi22y) + ux .* (E + p - tau11 - xi11) + vy .* (E + p - tau22 - xi22) - (tau12 + xi12) .* (vx + uy) + q1x + q2y + jE1x + jE2y];
   
end


