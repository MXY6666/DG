% Get one entry of the Jacobian of the flux function with respect to the
% conserved variables, w, r, z and q
function [ent, flag] = getJacobianEntry(rho, u, E, p, w, r, z, tau, xi, dp, i, j, tc)

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
            case 6
                ent = zeros(size(rho));
                flag = 0;
            case 7
                ent = zeros(size(rho));
                flag = 0;
            otherwise
                error('Wrong index')
        end
    case 2
        switch j
            case 1
                ent = dp(:, :, 1) - u.^2 - r / tc.We; 
                flag = 2;
            case 2
                ent = 2 * u + dp(:, :, 2);
                flag = 2;
            case 3
                ent = dp(:, :, 3);
                flag = 2;
            case 4
                ent = dp(:, :, 4) + w / tc.We;
                flag = 2;
            case 5
                ent = -1 / tc.We * rho; 
                flag = 2;
            case 6
                ent = repmat(-4 / 3 / tc.Re, size(rho));
                flag = 1;
            case 7
                ent = zeros(size(rho));
                flag = 0;                
            otherwise
                error('Wrong index')
        end        
    case 3 
        switch j
            case 1
                ent = (dp(:, :, 1) - r / tc.We - (E + p - tau - xi) ./ rho) .* u + w .* z ./ tc.We;
                flag = 2;
            case 2
                ent = dp(:, :, 2) .* u + (E + p - tau - xi) ./ rho;
                flag = 2;
            case 3
                ent = (1 + dp(:, :, 3)) .* u;
                flag = 2;
            case 4
                ent = (dp(:, :, 4) + w / tc.We) .* u + rho .* z / tc.We;
                flag = 2;
            case 5
                ent = -1 / tc.We * rho .* u;
                flag = 2;
            case 6
                ent = -4 / 3 / tc.Re * u + rho .* w / tc.We;
                flag = 2;
            case 7
                ent = ones(size(rho));
                flag = 1;
            otherwise
                error('Wrong index')
        end  
    otherwise
        error('Wrong index')
end

end