% fills the eight boundary data arrays based on the boundary 
% parameterizations for the given domain idmn
function [xb_bdry, yb_bdry, xt_bdry, yt_bdry, xl_bdry, yl_bdry, xr_bdry, yr_bdry]...
        = b_parm(m, n, idmn, a1, b1, a2, b2, a3, b3, a4, b4, tb, tt, tl, tr)
  
% 1 = Square 
% 2 = Rectangle 
% 3 = Parallelogram 
% 4 = Trapezoid
% 5 = Quadrilateral  
if (idmn <= 5)
    xb_bdry = tb * a2 + (1 - tb) * a1;
    yb_bdry = tb * b2 + (1 - tb) * b1;
    
    xt_bdry = tt * a3 + (1 - tt) * a4;
    yt_bdry = tt * b3 + (1 - tt) * b4;  
    
    xl_bdry = tl * a4 + (1 - tl) * a1;
    yl_bdry = tl * b4 + (1 - tl) * b1;
    
    xr_bdry = tr * a3 + (1 - tr) * a2;
    yr_bdry = tr * b3 + (1 - tr) * b2;   
end

% Annulus/Horseshoe
if (idmn == 6)
    xb_bdry = -a1 * cos(pi * tb);
    yb_bdry =  a2 * sin(pi * tb);
    
    xt_bdry = -a3 * cos(pi * tt);
    yt_bdry =  a4 * sin(pi * tt);

    xl_bdry = -(1 - tl) * a1 - tl * a3;
    yl_bdry = zeros(1, n + 1);
    
    xr_bdry = (1 - tr) * a1 + tr * a3;
    yr_bdry = zeros(1, n + 1);  
end

% Swan 
if (idmn == 7) 
    xb_bdry = tb;
    yb_bdry = zeros(1, m + 1);
    
    xt_bdry = tt;
    yt_bdry = 1 - 4 * (1 - a2) * tt + 4 * (1 - a2) .* tt .* tt;

    xl_bdry = zeros(1, n + 1);
    yl_bdry = tl;
    
    xr_bdry = 1 + 4 * a1 * tr - 4 * a1 * tr .* tr;
    yr_bdry = tr;
end 

% Dome/Valley
if (idmn == 8)
    xb_bdry = tb;
    yb_bdry = zeros(1, m + 1);
    
    xt_bdry = tt;
    yt_bdry = 1 + 16 * (a1 - 1) * tt .* tt .* (1 - tt) .* (1 - tt);

    xl_bdry = zeros(1, n + 1);
    yl_bdry = tl;
    
    xr_bdry = ones(1, n + 1);
    yr_bdry = tr;
end 

% Airfoil
if (idmn == 9) 
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    ind = find(tb <= 0.1);
    xb_bdry(ind) = 10 * tb(ind);
    yb_bdry(ind) = 1;
    ind = find(tb > 0.1 & tb <= 0.5);
    q = (tb(ind) - 0.1) .* (tb(ind) - 0.2) .* (tb(ind) - 0.5);
    xb_bdry(ind) = -1 / 3 + 15 * tb(ind) - 50 / 3 * tb(ind) .* tb(ind) + 13.889 * q;
    yb_bdry(ind) = (1 - 5 * a1 / 3) + 20 * a1 * tb(ind) - 100 * a1 / 3 * tb(ind) .* tb(ind) + 25 * (40 * a1 - 3) / 9 * q;
    if (a1 == 0)
        sss = 1.25 * (tb(ind) - 0.1);
        xb_bdry(ind) = 2 - cos(2 * pi * sss);
        yb_bdry(ind) = 1 + 0.025 * sin(2 * pi * sss);
    end
    ind = find(tb > 0.5 & tb < 0.9);
    s = 1 - tb(ind);
    q = (s - 0.1) .* (s - 0.2) .* (s - 0.5);
    xb_bdry(ind) = -1 / 3 + 15 * s - 50 / 3 * s .* s + 13.889 * q;
    yb_bdry(ind) = (1 - 5 * a1 / 3) + 20 * a1 * s - 100 * a1 / 3 * s .* s + 25 * (40 * a1 - 3) / 9 * q;
    yb_bdry(ind) = 2 - yb_bdry(ind);
    if (a1 == 0)
        sss = 1.25 * (tb(ind) - 0.1);
        xb_bdry(ind) = 2 - cos(2 * pi * sss);
        yb_bdry(ind) = 1 + 0.025 * sin(2 * pi * sss);        
    end
    ind = find(tb >= 0.9);
    xb_bdry(ind) = 10 * (1 - tb(ind));
    yb_bdry(ind) = 1;

    xt_bdry = zeros(1, m + 1);
    yt_bdry = zeros(1, m + 1);
    ind = find(tt <= 0.4);
    xt_bdry(ind) = 10 * tt(ind);
    yt_bdry(ind) = 2;
    ind = find(tt > 0.4 & tt < 0.6);
    xt_bdry(ind) = 4;
    yt_bdry(ind) = 6 - 10 * tt(ind);
    ind = find(tt >= 0.6);
    xt_bdry(ind) = 10 * (1 - tt(ind));
    yt_bdry(ind) = 0;
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = 1 + tl;
    
    xr_bdry = zeros(1, n + 1);
    yr_bdry = 1 - tr;      
end

% Chevron
if (idmn == 10) 
    xb_bdry = tb;
    yb_bdry = zeros(1, m + 1);
    ind = find(tb <= 0.5);
    yb_bdry(ind) = 2 * a1 * tb(ind);
    ind = find(tb > 0.5);
    yb_bdry(ind) = 2 * a1 * (1 - tb(ind));
    
    xt_bdry = tt;
    yt_bdry = zeros(1, m + 1);
    ind = find(tt <= 0.5);
    yt_bdry(ind) = 1 + 2 * a2 * tt(ind);
    ind = find(tt > 0.5);
    yt_bdry(ind) = 1 + 2 * a2 * (1 - tt(ind));
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = tl;
    
    xr_bdry = ones(1, n + 1);
    yr_bdry = tr;
end

% Horn
if (idmn == 11) 
    xb_bdry = 3 * tb;
    yb_bdry = zeros(1, m + 1);
    yb_bdry(tb <= 2 / 3) = -1;
    ind = find(tb > 2 / 3);
    yb_bdry(ind) = 1 - 3 * tb(ind);
    
    xt_bdry = 3 * tt;
    yt_bdry = zeros(1, m + 1);
    yt_bdry(tt <= 2 / 3) = 1;
    ind = find(tt > 2 / 3);
    yt_bdry(ind) = 3 * tt(ind) - 1;
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = 2 * tl - 1;
    
    xr_bdry = 3 + 4 * (a1 - 3) * tr - 4 * (a1 - 3) * tr .* tr;
    yr_bdry = -2 + 4 * tr;
end

% Backstep
if (idmn == 12) 
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    ind = find(tb <= a1 / (4 + b1));
    xb_bdry(ind) = (4 + b1) * tb(ind);
    yb_bdry(ind) = b1;
    ind = find(tb > a1 / (4 + b1) & tb < (a1 + b1) / (4 + b1));
    xb_bdry(ind) = a1;
    yb_bdry(ind) = a1 + b1 - (4 + b1) * tb(ind); 
    ind = find(tb >= (a1 + b1) / (4 + b1));
    xb_bdry(ind) = (4 + b1) * tb(ind) - b1;
    
    xt_bdry = 4 * tt;
    yt_bdry = 2 * ones(1, m + 1);
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = (2 - b1) * tl + b1;
    
    xr_bdry = 4 * ones(1, n + 1);
    yr_bdry = 2 * tr;
end

% Plow
if (idmn == 13) 
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    ind = find(tb <= a1 / (2 + a1));
    xb_bdry(ind) = 2;
    yb_bdry(ind) = (2 + a1) * tb(ind);
    ind = find(tb > a1 / (2 + a1));
    xb_bdry(ind) = (2 + a1) * tb(ind) + (2 - a1);
    yb_bdry(ind) = a1;
    
    xt_bdry = zeros(1, m + 1);
    yt_bdry = zeros(1, m + 1);
    ind = find(tt <= b1);
    xt_bdry(ind) = 2 - 2 * cos(pi * tt(ind) / b1 / 2);
    yt_bdry(ind) = 2 * sin(pi * tt(ind) / b1 / 2);
    ind = find(tt > b1);
    xt_bdry(ind) = 2 * tt(ind) / (1 - b1) + 2 * (1 - 2 * b1) / (1 - b1);
    yt_bdry(ind) = 2;
    
    xl_bdry = 2 * (1 - tl);
    yl_bdry = zeros(1, n + 1);
    
    xr_bdry = 4 * ones(1, n + 1);
    yr_bdry = a1 + (2 - a1) * tr;    
end

% C-grid
if (idmn == 14)
    xlb = 2 * (a1 + b1);
    xlt = 2 * (a2 + b2);
    rb1 = a1  / xlb;
    rb2 = (a1 + 2 * b1) / xlb;
    rt1 = a2 / xlt;
    rt2 = (a2 + 2 * b2) / xlt;
    
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);    
    ind = find(tb <= rb1);
    xb_bdry(ind) = xlb * tb(ind);
    yb_bdry(ind) = b1;  
    ind = find(tb > rb1 & tb < rb2);
    xb_bdry(ind) = a1;
    yb_bdry(ind) = xlb * (0.5 - tb(ind));  
    ind = find(tb >= rb2);
    xb_bdry(ind) = xlb * (1 - tb(ind));
    yb_bdry(ind) = -b1;
    
    xt_bdry = zeros(1, m + 1);
    yt_bdry = zeros(1, m + 1);
    ind = find(tt <= rt1);
    xt_bdry(ind) = xlt * tt(ind);
    yt_bdry(ind) = b2;   
    ind = find(tt > rt1 & tt < rt2);
    xt_bdry(ind) = a2;
    yt_bdry(ind) = xlt * (0.5 - tt(ind));  
    ind = find(tt >= rt2);
    xt_bdry(ind) = xlt * (1 - tt(ind));
    yt_bdry(ind) = -b2;   
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = (1 - tl) * b1 + tl * b2;
    
    xr_bdry = zeros(1, n + 1);
    yr_bdry = -yl_bdry;   
end

% New Airfoil
if (idmn == 15) 
    xlt = 2 * (a2 + b2);
    t1t = a2 / xlt;
    t2t = (a2 + 2 * b2) / xlt;
    xlb = 2 * (a3 + 2 * b3);
    t1b = a1 / xlb;
    t2b = a3 / xlb;
    t3b = 0.5;
    t4b = 1 - t2b;
    t5b = 1 - t1b;
    r1 = (a3 - b3 - a1) / (t2b - t1b);
    r2 = (a3 - b3 - a1) / (t4b - t5b);
    r3 = -b3 / (a3 - b3 - a1) / (a3 - b3 - a1);  
    
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);    
    ind = find(tb <= t1b);
    xb_bdry(ind) = a1 / t1b * tb(ind);
    yb_bdry(ind) = 0;
    ind = find(tb > t1b & tb <= t2b);
    xb_bdry(ind) = r1 * tb(ind) + a1 - t1b * r1;
    yb_bdry(ind) = r3 * (xb_bdry(ind) - a1) .* (xb_bdry(ind) - 2 * a3 + 2 * b3 + a1);    
    ind = find(tb > t2b & tb <= t3b);
    arg = (pi / 2) * (t3b - tb(ind)) / (t3b - t2b);
    xb_bdry(ind) = (a3 - b3) + b3 * cos(arg);
    yb_bdry(ind) = b3 * sin(arg);    
    ind = find(tb > t3b & tb <= t4b);
    arg = (pi / 2) * (t3b - tb(ind)) / (t3b - t4b);
    xb_bdry(ind) = (a3 - b3) + b3 * cos(arg);
    yb_bdry(ind) = -b3 * sin(arg);    
    ind = find(tb > t4b & tb <= t5b);
    xb_bdry(ind) = r2 * tb(ind) + a1 - t5b * r2;
    yb_bdry(ind) = -r3 * (xb_bdry(ind) - a1) .* (xb_bdry(ind) - 2 * a3 + 2 * b3 + a1);    
    ind = find(tb > t5b);
    xb_bdry(ind) = a1 * (1 - tb(ind)) / (1 - t5b);
    yb_bdry(ind) = 0;  
    
    xt_bdry = zeros(1, m + 1);
    yt_bdry = zeros(1, m + 1);
    ind = find(tt < t1t);
    xt_bdry(ind) = xlt * tt(ind);
    yt_bdry(ind) = b2;
    ind = find(tt >= t1t & tt <= t2t);
    xt_bdry(ind) = a2;
    yt_bdry(ind) = xlt * (0.5 - tt(ind));
    ind = find(tt > t2t);
    xt_bdry(ind) = xlt * (1 - tt(ind));
    yt_bdry(ind) = -b2;
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = b2 * tl;
    xr_bdry = zeros(1, n + 1);
    yr_bdry = -yl_bdry;
end

% Tie
if (idmn == 16)
    a4 = a2 + a3 - a1;
 
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    ind = find(tb <= 0.5);
    xb_bdry(ind) = 2 * a2 * tb(ind);
    yb_bdry(ind) = a1 * (1 - 4 * tb(ind) .* tb(ind));
    ind = find(tb > 0.5);
    xb_bdry(ind) = 2 * a3 * tb(ind) + a2 -a3;
    yb_bdry(ind) = 0;
    
    xt_bdry = a4 * tt;
    yt_bdry = (a1 + a4) * ones(1, m + 1);
    
    xl_bdry = zeros(1, n + 1);
    yl_bdry = a4 * tl + a1;
    
    xr_bdry = zeros(1, n + 1);
    yr_bdry = zeros(1, n + 1);
    ind = find(tr <= 0.5);
    xr_bdry(ind) = a2 + a3;
    yr_bdry(ind) = 2 * a3 * tr(ind);
    ind = find(tr > 0.5);
    xr_bdry(ind) = a4 + 4 * a1 * (1 - tr(ind)) .* (1 - tr(ind));
    yr_bdry(ind) = 2 * (a1 + a4 - a3) * tr(ind) - (a1 + a4 - 2 * a3);
end

% Shell
if (idmn == 17) 
    a3 = a1 + a2;
    
    xb_bdry = a1 * (2 * tb - 1);
    yb_bdry = zeros(1, m + 1);
    
    xt_bdry = -a3 * cos(pi * tt);
    yt_bdry = a2 + a3 * sin(pi * tt);
    
    xl_bdry = -a3 + a2 * cos(0.5 * pi * tl);
    yl_bdry = a2 * sin(0.5 * pi * tl);
    
    xr_bdry = a3 + a2 * cos(pi * (1 - tr / 2));
    yr_bdry = a2 * sin(pi * (1 - tr / 2));
end

% New Horseshoe
if (idmn == 18) 
    xb_bdry = -cos(pi * tb);
    yb_bdry = sin(pi * tb);
    
    xt_bdry = -2 * cos(pi * tt);
    yt_bdry = 2 * a1 * sin(pi * tt);
    
    xl_bdry = -(1 + tl);
    yl_bdry = zeros(1, n + 1);
    
    xr_bdry = 1 + tr;
    yr_bdry = zeros(1, n + 1);
end

% Ring
if (idmn == 19)
    xb_bdry = a1 * cos(2 * pi * (1 - tb));
    yb_bdry = a1 * sin(2 * pi * (1 - tb));
    
    xt_bdry = a2 * cos(2 * pi * (1 - tt));
    yt_bdry = a2 * sin(2 * pi * (1 - tt));    
    
    xl_bdry = a1 + (a2 - a1) * tl;
    yl_bdry = zeros(1, n + 1);
    
    xr_bdry = a1 + (a2 - a1) * tr;
    yr_bdry = zeros(1, n + 1);
end
    
% NACA0012 O-type
if (idmn == 20)
    % rf https://turbmodels.larc.nasa.gov/naca0012numerics_val.html
    coe = 0.594689181 * [0.298222773, -0.127125232, -0.357907906, 0.291984971, -0.105174606];
    
    frac = 0.5;
    ind = find(tb <= frac, 1, 'last');
    
    xb_bdry = 1 - cos(pi * (tb - 0.5));
    yb_bdry = zeros(1, m + 1);
    yb_bdry(1 : ind) = -(coe(1) * sqrt(xb_bdry(1 : ind)) + polyval([coe(5 : -1 : 2), 0], xb_bdry(1 : ind)));
    yb_bdry(ind + 1 : end) = coe(1) * sqrt(xb_bdry(ind + 1 : end)) + polyval([coe(5 : -1 : 2), 0], xb_bdry(ind + 1 : end));
    
    xt_bdry = a1 * cos(2 * pi * (1 - tt)) + 0.5;
    yt_bdry = a1 * sin(2 * pi * (1 - tt));   
    
    xl_bdry = 1 + (a1 - 0.5) * tl;
    yl_bdry = zeros(1, n + 1);
    
    xr_bdry = 1 + (a1 - 0.5) * tr;
    yr_bdry = zeros(1, n + 1);    
end

% NACA0012 C-type
if (idmn == 21)
    % rf https://turbmodels.larc.nasa.gov/naca0012numerics_val.html
    coe = 0.594689181 * [0.298222773, -0.127125232, -0.357907906, 0.291984971, -0.105174606];
       
    frac1 = (1 - b1) / 2;
    frac2 = 0.5;
    frac3 = (1 + b1) / 2;
    ind1 = find(tb >= frac1, 1);
    if (abs(tb(ind1 - 1) - frac1) < abs(tb(ind1) - frac1))
        ind1 = ind1 - 1;
    end
    tb(ind1) = frac1;
    ind2 = find(tb <= frac2, 1, 'last');
    ind3 = find(tb >= frac3, 1);
    if (abs(tb(ind3 - 1) - frac3) < abs(tb(ind3) - frac3))
        ind3 = ind3 - 1;
    end
    tb(ind3) = frac3;
    
    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    xb_bdry(1 : ind1) = a2 + (1 - a2) * tanh(2 * linspace(0, 1, ind1)) / tanh(2);
    xb_bdry(ind1 + 1 : ind3) = 1 - cos(pi / b1 * (tb(ind1 + 1 : ind3) - 0.5));
    xb_bdry(ind3 + 1 : end) = xb_bdry(ind1 - 1 : -1 : 1);
    yb_bdry(ind1 + 1 : ind2) = -(coe(1) * sqrt(xb_bdry(ind1 + 1 : ind2)) + polyval([coe(5 : -1 : 2), 0], xb_bdry(ind1 + 1 : ind2)));
    yb_bdry(ind2 + 1 : ind3) = coe(1) * sqrt(xb_bdry(ind2 + 1 : ind3)) + polyval([coe(5 : -1 : 2), 0], xb_bdry(ind2 + 1 : ind3));
    yb_bdry(ind3) = 0;
    
    frac1 = a2 / (2 * a2 + pi * a1);
    frac2 = (a2 + pi * a1) / (2 * a2 + pi * a1);
    ind1 = find(tt <= frac1, 1, 'last');
    ind2 = find(tt <= frac2, 1, 'last');
    
    xt_bdry = zeros(1, m + 1);
    yt_bdry = zeros(1, m + 1);      
    xt_bdry(1 : ind1) = a2 - (2 * a2 + pi * a1) * tt(1 : ind1);
    xt_bdry(ind1 + 1 : ind2) = -a1 * cos((2 * a2 + pi * a1) / a1 * tt(ind1 + 1 : ind2) - a2 / a1 - pi / 2);
    xt_bdry(ind2 + 1 : end) = (2 * a2 + pi * a1) * tt(ind2 + 1 : end) - (a2 + pi * a1);
    yt_bdry(1 : ind1) = -a1;
    yt_bdry(ind1 + 1 : ind2) = a1 * sin((2 * a2 + pi * a1) / a1 * tt(ind1 + 1 : ind2) - a2 / a1 - pi / 2);
    yt_bdry(ind2 + 1 : end) = a1;
    
    xl_bdry = a2 * ones(1, n + 1);
    yl_bdry = -a1 * tl;
    
    xr_bdry = a2 * ones(1, n + 1);
    yr_bdry = a1 * tr;
end

% Slope
if (idmn == 22)
    frac = (a2 - a1) / (a3 - a1);
    ind = find(tb >= frac, 1);
    if (abs(tb(ind - 1) - frac) < abs(tb(ind) - frac))
        ind = ind - 1;
    end
    tb(ind) = frac;

    xb_bdry = zeros(1, m + 1);
    yb_bdry = zeros(1, m + 1);
    xb_bdry(1 : ind) = a1 + (a2 - a1) / tb(ind) * tb(1 : ind);
    xb_bdry(ind : end) = a2 + (a3 - a2) / (1 - tb(ind)) * (tb(ind : end) - tb(ind));
    yb_bdry(1 : ind) = b1;
    yb_bdry(ind : end) = b1 + (b3 - b1) / (1 - tb(ind)) * (tb(ind : end) - tb(ind));

    xt_bdry = a1 + (a3 - a1) * tt;
    yt_bdry = b4 * ones(1, length(tt));

    xl_bdry = a1 * ones(1, length(tl));
    yl_bdry = b1 + (b4 - b1) * tl;

    xr_bdry = a3 * ones(1, length(tr));
    yr_bdry = b3 + (b4 - b3) * tr;
end

end


