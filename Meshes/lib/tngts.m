function  [xxi, yxi, xet, yet] = tngts(x, y, i, j, tdx, tdy)
 
% computes components of local tangents using centered-
% difference approximation

  xxi = (x(i + 1, j) - x(i - 1, j) ) / tdx;
  xet = (x(i, j + 1) - x(i, j - 1) ) / tdy;
  yxi = (y(i + 1, j) - y(i - 1, j) ) / tdx;
  yet = (y(i, j + 1) - y(i, j - 1) ) / tdy;  
  
end