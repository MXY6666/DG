function [p11,p12,p22,p11p11,p11p12,p11p22,p11prg,p12p12,p12p22,p12prg,p22p22,p22prg,prgprg]...
          = partials(g11, g12, g22, rg, g, imth)
 
p11 = 0.;
p12 = 0.;
p22 = 0.;

p11p11 = 0.;
p11p12 = 0.;
p11p22 = 0.;
p11prg = 0.;
p12p12 = 0.;
p12p22 = 0.;
p12prg = 0.;
p22p22 = 0.;
p22prg = 0.;
prgprg = 0.;

% ao
if (imth == 1)  
    p11 = g22;
    p22 = g11;
    p11p22 = 1;
    return
end

% ttm
if (imth == 2)  
    p11 = 1 / rg;
    p22 = 1 / rg;
    p11prg = -1 / g;
    p22prg = -1 / g;
    prgprg = 2 * (g11 + g22) / g / rg;
    return
end
% length
if (imth == 3)  
    p11 = 1;
    p22 = 1;       
    return
end
% area
if (imth == 4)  
    p11 = g22;
    p12 = -2 * g12;
    p22 = g11;
    p11p22 = 1;
    p12p12 = -2;
    return
end
% 0.9 area + 0.1 length
if (imth == 5)  
    w = 0.1;
    p11 = w * 1;
    p22 = w * 1;
    prgprg = (1 - w) * 2;
    return
end
% g12*g12
if (imth == 6)  
    p12 = 2 * g12;
    p12p12 = 2;
    return
end
% Liao
if (imth == 7)  
    p11 = 2 * g11;
    p12 = 4 * g12;
    p22 = 2 * g22;
    p11p11 = 2;
    p12p12 = 4;
    p22p22 = 2;
    return
end
% mod Liao
if (imth == 8)  
    ca = (g11 + g22) / rg;
    p11 = 2 * ca / rg;
    p22 = p11;
    p11p11 = 2 / g;
    p11p22 = p11p11;
    p22p22 = p11p11;
    p11prg = -4 * ca / g;
    p22prg = p11prg;
    prgprg = 6 * ca * ca / g;
    return
end
% scaled-laplace
if (imth == 9)  
    r = sqrt(g22 / g11);
    p11 = r / 2;
    p22 = 0.5 / r;
    p11p11 = -0.25 * r / g11;
    p11p22 = 0.25 / sqrt(g11 * g22);
    p22p22 = -0.25 / r / g22;
    return
end
% g12*g12/g11/g22
if (imth == 10)  
    p11 = -g12 * g12 / g11 / g11 /g22;
    p12 = 2 * g12 / g11 / g22;
    p22 = -g12 * g12 / g11 / g22 / g22;
    p11p11 = 2 * g12 * g12 / g11 / g11 / g11 / g22;
    p11p12 = -2 * g12 / g11 / g11 / g22;
    p11p22 = g12 * g12 / g11 / g11 / g22 / g22;
    p12p12 = 2 / g11 / g22;
    p12p22 = -2 * g12 / g11 / g22 / g22;
    p22p22 = 2 * g12 * g12 / g11 / g22 / g22 / g22;
    return
end
% g11+g22 - rg*(r + 1/r)
if (imth == 11)  
    r = sqrt(g22 / g11);
    s = 1 / r;
    u = sqrt(g / g11 / g22);
    p11 = 1 - u + 0.5 * u * (1 + r * r);
    p22 = 1 - u + 0.5 * u * (1 + s * s);
    p11p11 = (g11 - 3 * g22) * u / 4 / g11 / g11;
    p11p22 = 0.25 * u * (1 / g11 + 1 / g22);
    p22p22 = (g22 - 3 * g11) * u / 4 / g22 / g22;
    p11prg = (r - s) / 2 / g11;
    p22prg = (s - r) / 2 / g22;
    return
end
% AO**2
if (imth == 12)  
    p11 = 2 * g11 * g22 * g22;
    p22 = 2 * g22 * g11 * g11;
    p11p11 = 2 * g22 * g22;
    p11p22 = 4 * g11 * g22;
    p22p22 = 2 * g11 * g11;
end
% Guest: g + g12*g12/g11/g22  (modified AO)
if (imth == 13)  
%
%         p11 = g22 - g12*g12/g11/g11/g22
%         p12 = 2.*g12*(1./g11/g22 -1.)
%         p22 = g11 - g12*g12/g22/g22/g11
%         p11p11 = 2.*g12*g12/g11/g11/g11/g22
%         p11p12 = -2.*g12/g11/g11/g22
%         p11p22 = 1. + g12*g12/g11/g11/g22/g22
%         p12p12 = 2.*(1./g11/g22 - 1.)
%         p22p22 = 2.*g12*g12/g22/g22/g22/g11
%         p12p22 = -2.*g12/g22/g22/g11
%
%        Guest:  ( sqrt( g11 ) + sqrt( g22 ) )^2
%
p11 = 1 + sqrt(g22 / g11);
p22 = 1 + sqrt(g11 / g22 );
p11p11 = -0.5 * sqrt(g22 / g11) / g11;
p11p22 = 0.5 / sqrt(g11 * g22);
p22p22 = -0.5 * sqrt(g11 / g22) / g22;

end

end