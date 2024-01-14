% the maximum absoulute eigenvalue 
function alpha = getEigenMax(rho, u, tc)

dp = rho .* getFreeEnergyDerivative2(rho, tc);

alpha = max(abs(u - sqrt(dp)), abs(u + sqrt(dp)));

end