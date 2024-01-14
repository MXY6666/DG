% Get one entry of the Jacobian of the flux function
function ent = getJacobianEntry(rho, u, E, i, j, tc)

switch i
    case 1
        switch j
            case 1
                ent = zeros(size(rho));
            case 2
                ent = ones(size(rho));
            case 3
                ent = zeros(size(rho));
            otherwise
                error('Wrong index')
        end
    case 2
        switch j
            case 1
                ent = 0.5 * (tc.gamma - 3) * u.^2; 
            case 2
                ent = (3 - tc.gamma) * u;
            case 3
                ent = repmat(tc.gamma - 1, size(rho));
            otherwise
                error('Wrong index')
        end        
    case 3 
        switch j
            case 1
                ent = (tc.gamma - 1) * u.^3 - tc.gamma * u .* E ./ rho;
            case 2
                ent = tc.gamma * E ./ rho - 1.5 * (tc.gamma - 1) * u.^2;
            case 3
                ent = tc.gamma * u;
            otherwise
                error('Wrong index')
        end  
    otherwise
        error('Wrong index')
end

end