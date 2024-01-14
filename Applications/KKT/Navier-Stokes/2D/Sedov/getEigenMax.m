function alpha = getEigenMax(rho, un, p, tc)

alpha = abs(un) + getSpeedOfSound(rho, p, tc);

end