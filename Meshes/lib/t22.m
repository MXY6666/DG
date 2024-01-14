function [c11, c12, c22] = t22(xm11, xm12, xm22, ...
          g11, g12, g22, rg, g,...
          p11, p12, p22, p11p11, p11p12, p11p22, p11prg,...
          p12p12, p12p22, p12prg, p22p22, p22prg, prgprg)
 
s11 = p12p12 / 2 - p12prg * g12 / rg + prgprg * g12 * g12 / g / 2;
s12 = p12p22 - p22prg * g12 / rg + p12prg * g11 / rg / 2 - prgprg * g12 * g11 / g / 2;
s22 = 2 * p22p22 + 2 * p22prg * g11 / rg + prgprg * g11 * g11 / g / 2;

c11 = p22 + s11 * xm11(1, 1) + s12 * xm12(1, 1) + s22 * xm22(1, 1);
c12 = s11 * xm11(1, 2) + s12 * xm12(1, 2) + s22 * xm22(1, 2);
c22 = p22 + s11 * xm11(2, 2) + s12 * xm12(2, 2) + s22 * xm22(2, 2);

end
