function [g11, g12, g22, rg, g] = mets(xxi, yxi, xet, yet)


g11 = xxi * xxi + yxi * yxi;
g12 = xxi * xet + yxi * yet;
g22 = xet * xet + yet * yet;
rg = xxi * yet - xet * yxi;
g = rg * rg;
      
end
