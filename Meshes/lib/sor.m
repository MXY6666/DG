%  iterates on a 9-point stencil in two independent variables, which
%  results from a pair of equations of the form :
%  a11*x11 + b11*x12 + c11*x22 + a12*y11 + b12*y12 + c12*y22 = rhsx
%  a12*x11 + b12*x12 + c12*x22 + a22*y11 + b22*y12 + c22*y22 = rhsy
function [x, y, xerr] = sor(x, y, rhsx, rhsy, dxx, dxy, dyy, tdx, tdy, m, n, ncntmax, omg, imth)

xerr = 0;
if (imth > 13)
    for ncnt = 1 : ncntmax
        for i = 2 : m
            for j = 2 : n
                [xxi,yxi,xet,yet] = tngts(x,y,i,j,tdx,tdy);
                J2 = (xxi * yet - xet * yxi)^2;
                
                alpha = xet^2 + yet^2;
                beta = -(xxi * xet + yxi * yet);
                gamma = xxi^2 + yxi^2;
                phi = -J2 * (rhsx(i, j) * xxi + rhsy(i, j) * xet);
                psi = -J2 * (rhsx(i, j) * yxi + rhsy(i, j) * yet);
                
                c22 = -2 * (alpha / dxx + gamma / dyy);
                c32 = alpha / dxx;
                c12 = c32;
                c23 = gamma / dyy;
                c21 = c23;
                c33 = beta / (2 * dxy);
                c13 = -c33;
                c11 = c33;
                c31 = -c33;
                
                r = c22 * x(i,j) + c32 * x(i+1,j) + c33 * x(i+1,j+1) ...
                  + c23 * x(i,j+1) + c13 * x(i-1,j+1) + c12 * x(i-1,j) ...
                  + c11 * x(i-1,j-1) + c21 * x(i,j-1) + c31 * x(i+1,j-1) - phi;
                tmp = x(i,j) - omg * r / c22;
                xerr = max(xerr, abs(tmp - x(i,j)));
                x(i,j) = tmp;
                
                r = c22 * y(i,j) + c32 * y(i+1,j) + c33 * y(i+1,j+1) ...
                  + c23 * y(i,j+1) + c13 * y(i-1,j+1) + c12 * y(i-1,j) ...
                  + c11 * y(i-1,j-1) + c21 * y(i,j-1) + c31 * y(i+1,j-1) - psi;
                tmp = y(i,j) - omg * r / c22;    
                xerr = max(xerr, abs(tmp - y(i,j)));
                y(i,j) = tmp;              
            end         
        end
    end 
    return
end

for ncnt = 1 : ncntmax
    for j = 2 : n
        for i = 2 : m
            [xxi,yxi,xet,yet] = tngts(x,y,i,j,tdx,tdy);
            [g11,g12,g22,rg,g] = mets(xxi,yxi,xet,yet);
            [xm11,xm12,xm22] = tnsr(xxi,yxi,xet,yet);

            [p11,p12,p22,p11p11,p11p12,p11p22,p11prg,...
            p12p12,p12p22,p12prg,p22p22,p22prg,prgprg]...
            = partials(g11,g12,g22,rg,g,imth);

            [a11,a12,a22] = t11(xm11,xm12,xm22,g11,g12,g22,rg,g,...
            p11,p12,p22,p11p11,p11p12,p11p22,p11prg,...
            p12p12,p12p22,p12prg,p22p22,p22prg,prgprg );

            [b11,b12,b22] = t12(xm11,xm12,xm22,g11,g12,g22,rg,g,...
            p11,p12,p22,p11p11,p11p12,p11p22,p11prg,...
            p12p12,p12p22,p12prg,p22p22,p22prg,prgprg);

            [c11,c12,c22] = t22(xm11,xm12,xm22,g11,g12,g22,rg,g,...
            p11,p12,p22,p11p11,p11p12,p11p22,p11prg,...
            p12p12,p12p22,p12prg,p22p22,p22prg,prgprg);

            cf11 = 2 * (a11 / dxx + c11 / dyy );
            cf12 = 2 * (a12 / dxx + c12 / dyy );
            cf22 = 2 * (a22 / dxx + c22 / dyy );

            x12 = (x(i+1,j+1) - x(i+1,j-1) - x(i-1,j+1) + x(i-1,j-1)) / dxy;
            y12 = (y(i+1,j+1) - y(i+1,j-1) - y(i-1,j+1) + y(i-1,j-1)) / dxy;

            x11 = (x(i+1,j) + x(i-1,j)) / dxx;
            y11 = (y(i+1,j) + y(i-1,j)) / dxx;
            x22 = (x(i,j+1) + x(i,j-1)) / dyy;
            y22 = (y(i,j+1) + y(i,j-1)) / dyy;

            rhs1 = a11*x11 + a12*y11 + b11*x12 + b12*y12 + c11*x22 + c12*y22 - rhsx(i,j);
            rhs2 = a12*x11 + a22*y11 + b12*x12 + b22*y12 + c12*x22 + c22*y22 - rhsy(i,j);

            den = cf11 * cf22 - cf12 * cf12;
            if  (den == 0) 
                break;
            end
            xtmp = (cf22 * rhs1 - cf12 * rhs2) / den;
            ytmp = (cf11 * rhs2 - cf12 * rhs1) / den;
            xtmp = x(i,j) + omg * (xtmp - x(i,j));
            ytmp = y(i,j) + omg * (ytmp - y(i,j));
            xerr = max(xerr, abs(xtmp - x(i,j)));
            xerr = max(xerr, abs(ytmp - y(i,j)));
            x(i,j) = xtmp;
            y(i,j) = ytmp;
        end
    end
end

end
