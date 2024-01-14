% speed of sound
function c = getSpeedOfSound(rho, p, tc)

c = sqrt(abs(tc.gamma * p ./ rho));

end