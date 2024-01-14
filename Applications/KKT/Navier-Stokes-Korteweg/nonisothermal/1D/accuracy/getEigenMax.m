% the maximum absoulute eigenvalue 
function alpha = getEigenMax(rho, u, T, p, tc)

beta = 2 * (p ./ rho + 4 * T .* (3 + tc.Cv * (2 * rho - 3)) ./ (tc.Cv * (rho - 3).^2));

alpha = max(abs(u - sqrt(beta)), abs(u + sqrt(beta)));

end