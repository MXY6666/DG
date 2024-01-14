% Get one entry of the Jacobian of the flux function in y direction with 
% respect to the conserved variables, w, r, z and q 
function [ent, flag] = getJacobianEntry2(rho, u, v, E, p, w1, w2, r, trz, tau12, tau22, xi12, xi22, dp, i, j, tc)

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
            case 11
                ent = zeros(size(rho));
                flag = 0;
            case 12
                ent = zeros(size(rho));
                flag = 0;
            case 13
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
                ent = w2 / tc.We; 
                flag = 2;
            case 6
                ent = w1 / tc.We; 
                flag = 2;
            case 7
                ent = zeros(size(rho)); 
                flag = 0;
            case 8
                ent = zeros(size(rho)); 
                flag = 0;
            case 9
                ent = repmat(-1 / tc.Re, size(rho));
                flag = 1;
            case 10
                ent = repmat(-1 / tc.Re, size(rho));
                flag = 1;
            case 11
                ent = zeros(size(rho)); 
                flag = 0;
            case 12
                ent = zeros(size(rho)); 
                flag = 0;
            case 13
                ent = zeros(size(rho)); 
                flag = 0;
            otherwise
                error('Wrong index')
        end           
    case 3
        switch j
            case 1
                ent = dp(:, :, 1) - v.^2 - r / tc.We;  
                flag = 2;
            case 2
                ent = dp(:, :, 2);
                flag = 2;
            case 3
                ent = 2 * v + dp(:, :, 3);
                flag = 2;
            case 4
                ent = dp(:, :, 4);
                flag = 2;
            case 5
                ent = dp(:, :, 5) - w1 / tc.We;
                flag = 2;
            case 6
                ent = dp(:, :, 6) + w2 / tc.We;
                flag = 2;
            case 7
                ent = -1 / tc.We * rho; 
                flag = 2;
            case 8
                ent = repmat(2 / 3 / tc.Re, size(rho));
                flag = 1;
            case 9
                ent = zeros(size(rho));    
                flag = 0;
            case 10
                ent = zeros(size(rho));    
                flag = 0;
            case 11
                ent = repmat(-4 / 3 / tc.Re, size(rho));
                flag = 1;
            case 12
                ent = zeros(size(rho));    
                flag = 0;
            case 13
                ent = zeros(size(rho));    
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 4
        switch j
            case 1
                ent = (dp(:, :, 1) - r / tc.We - (E + p - tau22 - xi22) ./ rho) .* v + (tau12 + xi12) .* u ./ rho  + w2 .* trz ./ tc.We;
                flag = 2;
            case 2
                ent = dp(:, :, 2) .* v - (tau12 + xi12) ./ rho;
                flag = 2;
            case 3
                ent = dp(:, :, 3) .* v + (E + p - tau22 - xi22) ./ rho;
                flag = 2;
            case 4
                ent = (1 + dp(:, :, 4)) .* v;
                flag = 2;
            case 5
                ent = (dp(:, :, 5) - w1 / tc.We) .* v + w2 .* u / tc.We;
                flag = 2;
            case 6 
                ent = (dp(:, :, 6) + w2 / tc.We) .* v + (w1 .* u + rho .* trz) / tc.We;
                flag = 2;
            case 7
                ent = -1 / tc.We * rho .* v;
                flag = 2;
            case 8
                ent = 2 / 3 / tc.Re * v + rho .* w2 / tc.We;
                flag = 2;
            case 9
                ent = -1 / tc.Re * u;
                flag = 2;
            case 10
                ent = -1 / tc.Re * u;
                flag = 2;
            case 11
                ent = -4 / 3 / tc.Re * v + rho .* w2 / tc.We;
                flag = 2;
            case 12
                ent = zeros(size(rho));
                flag = 0;
            case 13       
                ent = ones(size(rho));
                flag = 1;
            otherwise
                error('Wrong index')
        end        
    otherwise
        error('Wrong index')
end

end
