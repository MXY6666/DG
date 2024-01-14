function [b11, b12, b22] = t12(xm11, xm12, xm22, ...
          g11, g12, g22, rg, g,...
          p11, p12, p22, p11p11, p11p12, p11p22, p11prg,...
          p12p12, p12p22, p12prg, p22p22, p22prg, prgprg)
 
s11 = 2 * p11p12 - 2 * p11prg * g12 / rg - prgprg * g12 * g22 / g;
s12 = 2 * p11p22 + p11prg * g11 / rg + p12p12 / 2 + p22prg * g22 / rg + prgprg * (g11 * g22 + g12 * g12) / g / 2;
s22 = 2 * p12p22 - 2 * p22prg * g12 / rg - prgprg * g11 * g12 / g;

b11 = p12 + rg * p12prg + s11 * xm11(1, 1) + s12 * xm12(1, 1) + s22 * xm22(1, 1);
b12 = s11 * xm11(1, 2) + s12 * xm12(1, 2) + s22 * xm22(1, 2);
b22 = p12 + rg * p12prg + s11 * xm11(2, 2) + s12 * xm12(2, 2) + s22 * xm22(2, 2);

end