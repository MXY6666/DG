% Computes the following grid "badness" measures:
%  (i) folding,
% (ii) angle non-uniformity
%(iii) area non-uniformity
function grd_met(x, y, m, n)

fold = false;
angmax = 0; angmin = 1.e+30; angmean = 0; angsq = 0;
kangmx = 0; iangmx = 0; kangmn = 0; iangmn = 0;
armax = 0; armin = 1.e+30; armean = 0; arsq = 0;
karmx = 0; karmn = 0;
den = 4 * m * n - 4;
rd = 180 / pi;

angbar = zeros(1, 4 * m * n); 
arbar = zeros(1, m * n);
for k = 1 : m * n
    j = ceil(k / m);
    i = k - (j - 1) * m;
    
    [lsw, asw, arsw] = sw(i, j, x, y);
    [lse, ase, arse] = se(i, j, x, y);
    [lnw, anw, arnw] = nw(i, j, x, y);
    [lne, ane, arne] = ne_(i, j, x, y);
    
    angbar(4 * k - 3) = rd * asw;
    angbar(4 * k - 2) = rd * ase;
    angbar(4 * k - 1) = rd * anw;
    angbar(4 * k) = rd * ane;
    arbar(k) = (arsw + arse + arnw + arne) / 4;  
         
    lfold = folding(lsw, lse, lnw, lne, k);
    fold = lfold || fold;
    
    [angmax,kangmx,iangmx,angmin,kangmn,iangmn,angmean,angsq] = angles(asw,ase,anw,ane,k,m,n,angmax,kangmx,iangmx,angmin,kangmn,iangmn,angmean,angsq);
    
    [armax,karmx,armin,karmn,armean,arsq] = areas(arsw,arse,arnw,arne,k,armax,karmx,armin,karmn,armean,arsq);
end
if ~fold
    disp('Grid is unfolded')
end

disp('Angles:')
stats(rd,den,angmean,angsq,angmax,kangmx,iangmx,angmin,kangmn,iangmn);

disp('Areas:')
stats(1,den+4,armean,arsq/4,armax/4,karmx,karmx,armin/4,karmn,karmn);
    
end

function [lfld, ang, cp] = sw(i, j, x, y)

a = [x(i+1,j) - x(i,j), y(i+1,j) - y(i,j)];
b = [x(i,j+1) - x(i,j), y(i,j+1) - y(i,j)];

[lfld, ang, cp] = xyz(a, b);

end

function [lfld, ang, cp] = se(i, j, x, y)

a = [x(i+1,j+1) - x(i+1,j), y(i+1,j+1) - y(i+1,j)];
b = [x(i,j) - x(i+1,j), y(i,j) - y(i+1,j)];

[lfld, ang, cp] = xyz(a, b);

end

function [lfld, ang, cp] = ne_(i, j, x, y)

a = [x(i,j+1) - x(i+1,j+1), y(i,j+1) - y(i+1,j+1)];
b = [x(i+1,j) - x(i+1,j+1), y(i+1,j) - y(i+1,j+1)];

[lfld, ang, cp] = xyz(a, b);

end

function [lfld, ang, cp] = nw(i, j, x, y)

a = [x(i,j) - x(i,j+1), y(i,j) - y(i,j+1)];
b = [x(i+1,j+1) - x(i,j+1), y(i+1,j+1) - y(i,j+1)];

[lfld, ang, cp] = xyz(a, b);

end

function [lfld, ang, cp] = xyz(a, b)

pi2 = pi / 2.;
lfld = false;

cp = cross(a, b);
if (cp <= 0) 
   lfld = true;
end

dp = dot(a, b);
if (dp > 0) 
    ang = abs(pi2 - atan(cp / dp)); 
elseif (dp < 0) 
    ang = abs(pi2 + atan(cp / dp));
else 
    if (cp <= 0) 
        ang = pi;
    else
        ang = 0;
    end       
end

end

function cp = cross(a, b)

cp = a(1) * b(2) - a(2) * b(1);

end

function dp = dot(a, b)

dp = a(:)' * b(:);

end

function lfold = folding(lsw, lse, lnw, lne, k)

lfold = false;
if (lsw || lse || lnw || lne)
    lfold = true;
    fprintf('Cell number %d is folded.\n', k)
end 

end

function [angmax,kangmx,iangmx,angmin,kangmn,iangmn,angmean,angsq] = angles(asw,ase,anw,ane,k,m,n,angmax,kangmx,iangmx,angmin,kangmn,iangmn,angmean,angsq)
    
if (k ~= 1)
    if (asw > angmax) 
        angmax = asw;
        iangmx = 1;
        kangmx = k;
    end 
    if (asw < angmin)
        angmin = asw;
        iangmn = 1;
        kangmn = k;
    end 
    angmean = angmean + asw;
    angsq = angsq + asw * asw;
end 

if (k ~= m)
    if (ase > angmax) 
        angmax = ase;
        iangmx = 2;
        kangmx = k;
    end 
    if (ase < angmin) 
        angmin = ase;
        iangmn = 2;
        kangmn = k;
    end 
    angmean = angmean + ase;
    angsq = angsq + ase*ase;
end

if (k ~= m * n) 
    if (ane > angmax) 
        angmax = ane;
        iangmx = 3;
        kangmx = k;
    end 
    if (ane < angmin) 
        angmin = ane;
        iangmn = 3;
        kangmn = k;
    end 
    angmean = angmean + ane;
    angsq = angsq + ane * ane;
end 

if (k ~= m*(n-1)+1) 
    if (anw > angmax)
        angmax = anw;
        iangmx = 4;
        kangmx = k;
    end 
    if (anw < angmin) 
        angmin = anw;
        iangmn = 4;
        kangmn = k;
    end 
    angmean = angmean + anw;
    angsq = angsq + anw * anw;
end
   
end

function [armax,karmx,armin,karmn,armean,arsq] = areas(arsw,arse,arnw,arne,k,armax,karmx,armin,karmn,armean,arsq)

area = arsw + arse + arnw + arne;

if (area > armax) 
    armax = area;
    karmx = k;
end 
if (area < armin) 
    armin = area;
    karmn = k;
end 

armean = armean + area;
arsq = arsq + area*area;

end

function stats(xunit,den,xmean,xsq,xmax,kmx,imx,xmin,kmn,imn)

xmean = xmean / den;
xdev = xsq / den - xmean * xmean;
xmax = xunit * xmax;
xmin = xunit * xmin;
xmean = xunit * xmean;
xdev = xunit * sqrt(abs(xdev));
if (xdev ~= 0) 
    xfar = abs(xmax - xmean) / xdev;
end
fprintf(' Mean: %f\n', xmean)
fprintf(' Std. Dev: %f\n', xdev)
fprintf(' Max: %f, %d, %d\n', xmax,kmx,imx)
fprintf(' Min: %f, %d, %d\n', xmin,kmn,imn)
end





