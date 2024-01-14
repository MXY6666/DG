% the maximum absoulute eigenvalue 
function alpha = getEigenMax(rho, un, T, p, tc)

beta = 2 * (p ./ rho + 4 * T .* (3 + tc.Cv * (2 * rho - 3)) ./ (tc.Cv * (rho - 3).^2));

alpha = max(abs(un - sqrt(beta)), abs(un + sqrt(beta)));

end