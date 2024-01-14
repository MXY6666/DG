function src = getSourceTerm(x, y, t, tc)

rho   = tc.rho(x, y, t);
rhox  = tc.rho_x(x, y, t);
rhoy  = tc.rho_y(x, y, t);
rhot  = tc.rho_t(x, y, t);
rhoxx = tc.rho_xx(x, y, t);
rhoyy = tc.rho_yy(x, y, t);

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

p   = tc.p(x, y, t);
px  = tc.p_x(x, y, t);
py  = tc.p_y(x, y, t);
pt  = tc.p_t(x, y, t);
pxx = tc.p_xx(x, y, t);
pyy = tc.p_yy(x, y, t);

m1     = rho .* u;
m2     = rho .* v;
u2     = u.^2;
uv     = u .* v;
v2     = v.^2;
q2     = u2 + v2;
E      = p  / (tc.gamma - 1) + 0.5 * rho  .* q2;
Ex     = px / (tc.gamma - 1) + 0.5 * rhox .* q2 + m1 .* ux + m2 .* vx;
Ey     = py / (tc.gamma - 1) + 0.5 * rhoy .* q2 + m1 .* uy + m2 .* vy;
Et     = pt / (tc.gamma - 1) + 0.5 * rhot .* q2 + m1 .* ut + m2 .* vt;
tau11  = tc.mu * (4 / 3 * ux - 2 / 3 * vy);
tau22  = tc.mu * (4 / 3 * vy - 2 / 3 * ux);
tau12  = tc.mu * (uy + vx);
tau11x = tc.mu * (4 / 3 * uxx - 2 / 3 * vxy);
tau22y = tc.mu * (4 / 3 * vyy - 2 / 3 * uxy);
tau12x = tc.mu * (uxy + vxx);
tau12y = tc.mu * (uyy + vxy);
q1x    = tc.gamma * tc.mu / (1 - tc.gamma) / tc.Pr * (pxx - (2 * rhox .* px + p .* rhoxx) ./ rho + 2 * p .* rhox .* rhox ./ rho ./ rho) ./ rho;
q2y    = tc.gamma * tc.mu / (1 - tc.gamma) / tc.Pr * (pyy - (2 * rhoy .* py + p .* rhoyy) ./ rho + 2 * p .* rhoy .* rhoy ./ rho ./ rho) ./ rho;

src = [rhot + rho .* (ux + vy) + rhox .* u + rhoy .* v; ...
       rhot .* u + rho .* ut + rhox .* u2 + rhoy .* uv + m1 .* (2 * ux + vy) + m2 .* uy + px - tau11x - tau12y; ...
       rhot .* v + rho .* vt + rhox .* uv + rhoy .* v2 + m2 .* (ux + 2 * vy) + m1 .* vx + py - tau12x - tau22y; ...
       Et + u .* (Ex + px - tau11x - tau12y) + v .* (Ey + py - tau12x - tau22y) + ux .* (E + p - tau11) + vy .* (E + p - tau22) - tau12 .* (vx + uy) + q1x + q2y];
   
end


