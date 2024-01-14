% Get one entry of the Jacobian of the flux function in y direction
function ent = getJacobianEntry2(rho, u, v, E, i, j, tc)

switch i 
    case 1
        switch j
            case 1
                ent = zeros(size(rho));
            case 2
                ent = zeros(size(rho));
            case 3
                ent = ones(size(rho));
            case 4
                ent = zeros(size(rho));                
            otherwise
                error('Wrong index')
        end        
    case 2
        switch j
            case 1
                ent = -u .* v;
            case 2
                ent = v;
            case 3
                ent = u;
            case 4
                ent = zeros(size(rho)); 
            otherwise
                error('Wrong index')
        end           
    case 3
        switch j
            case 1
                ent = 0.5 * ((tc.gamma - 1) * u.^2 + (tc.gamma - 3) * v.^2); 
            case 2
                ent = (1 - tc.gamma) * u;
            case 3
                ent = (3 - tc.gamma) * v;
            case 4
                ent = repmat(tc.gamma - 1, size(rho));
            otherwise
                error('Wrong index')
        end        
    case 4
        switch j
            case 1
                ent = ((tc.gamma - 1) * (u.^2 + v.^2) - tc.gamma * E ./ rho) .* v;
            case 2
                ent = (1 - tc.gamma) * u .* v;
            case 3
                ent = tc.gamma * E ./ rho - 0.5 * (tc.gamma - 1) * (u.^2 + 3 * v.^2);
            case 4
                ent = tc.gamma * v;
            otherwise
                error('Wrong index')
        end        
    otherwise
        error('Wrong index')
end


end