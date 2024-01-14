% Get one entry of the Jacobian of the flux function in y direction
function [ent, flag] = getJacobianEntry2(rho, u, v, E, tau12, tau22, i, j, tc)

switch i 
    case 1
        switch j
            case 1
                ent = zeros(size(rho));
                flag = 0;
            case 2
                ent = zeros(size(rho));
                flag = 0;
            case 3
                ent = ones(size(rho));
                flag = 1;
            case 4
                ent = zeros(size(rho));  
                flag = 0;
            case 5
                ent = zeros(size(rho));    
                flag = 0;
            case 6
                ent = zeros(size(rho));    
                flag = 0;
            case 7
                ent = zeros(size(rho));    
                flag = 0;
            case 8
                ent = zeros(size(rho));    
                flag = 0;
            case 9
                ent = zeros(size(rho));    
                flag = 0;
            case 10
                ent = zeros(size(rho));    
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 2
        switch j
            case 1
                ent = -u .* v;
                flag = 2;
            case 2
                ent = v;
                flag = 2;
            case 3
                ent = u;
                flag = 2;
            case 4
                ent = zeros(size(rho)); 
                flag = 0;
            case 5
                ent = zeros(size(rho)); 
                flag = 0;
            case 6
                ent = repmat(-tc.mu, size(rho));
                flag = 1;
            case 7
                ent = repmat(-tc.mu, size(rho));
                flag = 1;
            case 8
                ent = zeros(size(rho)); 
                flag = 0;
            case 9
                ent = zeros(size(rho)); 
                flag = 0;
            case 10
                ent = zeros(size(rho)); 
                flag = 0;
            otherwise
                error('Wrong index')
        end           
    case 3
        switch j
            case 1
                ent = 0.5 * ((tc.gamma - 1) * u.^2 + (tc.gamma - 3) * v.^2); 
                flag = 2;
            case 2
                ent = (1 - tc.gamma) * u;
                flag = 2;
            case 3
                ent = (3 - tc.gamma) * v;
                flag = 2;
            case 4
                ent = repmat(tc.gamma - 1, size(rho));
                flag = 1;
            case 5
                ent = repmat(2 / 3 * tc.mu, size(rho));
                flag = 1;
            case 6
                ent = zeros(size(rho));    
                flag = 0;
            case 7
                ent = zeros(size(rho));    
                flag = 0;
            case 8
                ent = repmat(-4 / 3 * tc.mu, size(rho));
                flag = 1;
            case 9
                ent = zeros(size(rho));    
                flag = 0;
            case 10
                ent = zeros(size(rho));    
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 4
        switch j
            case 1
                ent = ((tc.gamma - 1) * (u.^2 + v.^2) - tc.gamma * E ./ rho) .* v + (tau12 .* u + tau22 .* v) ./ rho;
                flag = 2;
            case 2
                ent = (1 - tc.gamma) * u .* v - tau12 ./ rho;
                flag = 2;
            case 3
                ent = (tc.gamma * E - tau22) ./ rho - 0.5 * (tc.gamma - 1) * (u.^2 + 3 * v.^2);
                flag = 2;
            case 4
                ent = tc.gamma * v;
                flag = 2;
            case 5
                ent = 2 / 3 * tc.mu * v;
                flag = 2;
            case 6
                ent = -tc.mu * u;
                flag = 2;
            case 7
                ent = -tc.mu * u;
                flag = 2;
            case 8
                ent = -4 / 3 * tc.mu * v;
                flag = 2;
            case 9
                ent = zeros(size(rho));
                flag = 0;
            case 10       
                ent = ones(size(rho));
                flag = 1;
            otherwise
                error('Wrong index')
        end        
    otherwise
        error('Wrong index')
end

end
