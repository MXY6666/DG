function [src1, src2] = getSourceTerm(x, t, tc)

rho    = tc.rho(x, t);
rhot   = tc.rho_t(x, t);
rhox   = tc.rho_x(x, t);
rhoxxx = tc.rho_xxx(x, t);

u  = tc.u(x, t);
ut = tc.u_t(x, t);
ux = tc.u_x(x, t);

d2F = getFreeEnergyDerivative2(rho, tc);

src1 = rhot + rhox .* u + rho .* ux;
src2 = ut + u .* ux + d2F .* rhox - tc.kappa .* rhoxxx;

end