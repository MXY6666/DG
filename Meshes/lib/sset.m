% set the source term of original pde
function [rhsx, rhsy] = sset(m, n, idmn, tb, tt, tl, tr)

rhsx = zeros(m + 1, n + 1);
rhsy = zeros(m + 1, n + 1);
if (idmn == 20)
    rhsy = repmat(-exp(-50*tl), [m + 1, 1]);
end
if (idmn == 21)
    rhsy = repmat(-exp(-50*tl), [m + 1, 1]);
end