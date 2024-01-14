% the maximum absoulute value(|u| + c)of the three real eigenvalues
% (u - c, u, u + c) of the Jacobian matrix of the flux function 
function alpha = getEigenMax(rho, u, p, tc)

alpha = abs(u) + getSpeedOfSound(rho, p, tc);

end