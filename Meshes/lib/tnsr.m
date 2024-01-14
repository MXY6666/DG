function [xm11, xm12, xm22] = tnsr(xxi, yxi, xet, yet)

xm11 = zeros(2, 2);
xm12 = zeros(2, 2);
xm22 = zeros(2, 2);

xm11(1, 1) = xxi * xxi;
xm11(1, 2) = xxi * yxi;
xm11(2, 2) = yxi * yxi;

xm12(1, 1) = 2 * xxi * xet;
xm12(1, 2) = xxi * yet + yxi * xet;
xm12(2, 2) = 2 * yxi * yet;

xm22(1, 1) = xet * xet;
xm22(1, 2) = xet * yet;
xm22(2, 2) = yet * yet;
 
end

