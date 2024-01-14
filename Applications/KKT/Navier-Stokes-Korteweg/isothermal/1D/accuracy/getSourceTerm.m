function [src1, src2] = getSourceTerm(x, t, tc)

rho    = tc.rho(x, t);
rhot   = tc.rho_t(x, t);
rhox   = tc.rho_x(x, t);
rhoxxx = tc.rho_xxx(x, t);

u   = tc.u(x, t);
ut  = tc.u_t(x, t);
ux  = tc.u_x(x, t);
uxx = tc.u_xx(x, t);

d2F = getFreeEnergyDerivative2(rho, tc);

src1 = rhot + rhox .* u + rho .* ux;
src2 = rhot .* u + rho .* ut + 2 * rho .* u .* ux + (u.^2 + rho .* d2F) .* rhox - 4 / 3 / tc.Re * uxx - 1 / tc.We * rho .* rhoxxx;

end