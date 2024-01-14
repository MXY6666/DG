function src = getSourceTerm(x, t, tc)

rho   = tc.rho(x, t);
rhot  = tc.rho_t(x, t);
rhox  = tc.rho_x(x, t);
rhoxx = tc.rho_xx(x, t);

u   = tc.u(x, t);
ut  = tc.u_t(x, t);
ux  = tc.u_x(x, t);
uxx = tc.u_xx(x, t);

p   = tc.p(x, t);
pt  = tc.p_t(x, t);
px  = tc.p_x(x, t);
pxx = tc.p_xx(x, t);

m    = rho .* u;
u2   = u.^2;
E    = p  / (tc.gamma - 1) + 0.5 * rho  .* u2;
Ex   = px / (tc.gamma - 1) + 0.5 * rhox .* u2 + m .* ux;
Et   = pt / (tc.gamma - 1) + 0.5 * rhot .* u2 + m .* ut;
tau  = 4 / 3 * tc.mu * ux;
taux = 4 / 3 * tc.mu * uxx;
qx   = tc.gamma * tc.mu / (1 - tc.gamma) / tc.Pr * (pxx - (2 * rhox .* px + p .* rhoxx) ./ rho + 2 * p .* rhox .* rhox ./ rho ./ rho) ./ rho;

src = [rhot + rhox .* u + rho .* ux; ...
       rhot .* u + rho .* ut + rhox .* u2 + 2 * m .* ux + px - taux; ...
       Et + ux .* (E + p - tau) + u .* (Ex + px - taux) + qx];

end