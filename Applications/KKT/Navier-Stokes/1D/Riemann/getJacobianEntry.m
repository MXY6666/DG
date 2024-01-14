% Get one entry of the Jacobian of the flux function with respect to the
% conserved variables, z and q
function [ent, flag] = getJacobianEntry(rho, u, E, tau, i, j, isViscous, tc)

if (~isViscous && j > 3)
    error('Wrong index')
end

switch i
    case 1
        switch j
            case 1
                ent = zeros(size(rho));
                flag = 0;
            case 2
                ent = ones(size(rho));
                flag = 1;
            case 3
                ent = zeros(size(rho));
                flag = 0;
            case 4
                ent = zeros(size(rho));
                flag = 0;                
            case 5
                ent = zeros(size(rho));
                flag = 0;                
            otherwise
                error('Wrong index')
        end
    case 2
        switch j
            case 1
                ent = 0.5 * (tc.gamma - 3) * u.^2; 
                flag = 2;
            case 2
                ent = (3 - tc.gamma) * u;
                flag = 2;
            case 3
                ent = repmat(tc.gamma - 1, size(rho));
                flag = 1;
            case 4
                ent = repmat(-4 / 3 * tc.mu, size(rho));
                flag = 1;
            case 5
                ent = zeros(size(rho));
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 3 
        switch j
            case 1
                ent = (tc.gamma - 1) * u.^3 - tc.gamma * u .* E ./ rho;
                if isViscous
                    ent = ent + tau .* u ./ rho;
                end
                flag = 2;
            case 2
                ent = tc.gamma * E ./ rho - 1.5 * (tc.gamma - 1) * u.^2;
                if isViscous
                    ent = ent - tau ./ rho;
                end
                flag = 2;
            case 3
                ent = tc.gamma * u;
                flag = 2;
            case 4
                ent = -4 / 3 * tc.mu * u;
                flag = 2;
            case 5
                ent = ones(size(rho));
                flag = 1;
            otherwise
                error('Wrong index')
        end  
    otherwise
        error('Wrong index')
end

end