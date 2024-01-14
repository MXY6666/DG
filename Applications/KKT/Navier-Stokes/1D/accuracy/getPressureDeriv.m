% Compute the derivative of pressure with respect to the conserved
% variables
function dp = getPressureDeriv(u, i, tc)

switch i
    case 1
        dp = 0.5 * (tc.gamma - 1) * u.^2;
    case 2
        dp = (1 - tc.gamma) * u;
    case 3
        dp = repmat(tc.gamma - 1, size(u));
    otherwise
        error('Wrong index')
end

end