function src = getSourceTerm(x, t, tc)

rho    = tc.rho(x, t);
rhox   = tc.rho_x(x, t);
rhot   = tc.rho_t(x, t);
rhoxx  = tc.rho_xx(x, t);
rhoxt  = tc.rho_xt(x, t);
rhoxxx = tc.rho_xxx(x, t);

u   = tc.u(x, t);
ux  = tc.u_x(x, t);
ut  = tc.u_t(x, t);
uxx = tc.u_xx(x, t);

T   = tc.T(x, t);
Tx  = tc.T_x(x, t);
Tt  = tc.T_t(x, t);
Txx = tc.T_xx(x, t);

m    = rho .* u;
u2   = u.^2;
tmp  = 8 / 3 * tc.Cv * T - 3 * rho + 0.5 * u.^2;
E    = rho  .* tmp + 0.5 / tc.We * rhox.^2;
Ex   = rhox .* tmp + rho .* (8 / 3 * tc.Cv * Tx - 3 * rhox + u .* ux) + rhox .* rhoxx / tc.We;
Et   = rhot .* tmp + rho .* (8 / 3 * tc.Cv * Tt - 3 * rhot + u .* ut) + rhox .* rhoxt / tc.We;
p    = 8 * rho .* T ./ (3 - rho) - 3 * rho.^2;
px   = 24 * rhox ./ (3 - rho).^2 .* T + 8 * rho ./ (3 - rho) .* Tx - 6 * rho .* rhox;
tau  = 4 / 3 / tc.Re * ux;
taux = 4 / 3 / tc.Re * uxx;
xi   = 1 / tc.We * (rho .* rhoxx - 0.5 * rhox.^2);
xix  = 1 / tc.We * (rho .* rhoxxx);
qx   = -(8 * tc.Cv) / (3 * tc.We * tc.Pr) * Txx;
jEx  = 1 / tc.We * (rhox.^2 .* ux + rho .* rhoxx .* ux + rho .* rhox .* uxx);

src = [rhot + rhox .* u + rho .* ux; ...
       rhot .* u + rho .* ut + rhox .* u2 + 2 * m .* ux + px - taux - xix; ...
       Et + ux .* (E + p - tau - xi) + u .* (Ex + px - taux - xix) + qx + jEx];

end