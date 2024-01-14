% Compute the Hessian matrix of pressure with respect to the conserved
% variables
function dp2 = getPressureHessian(rho, u, v, E, w1, w2, i, j, tc)

switch i
    case 1
        switch j
            case 1
                dp2 = 6 / tc.Cv ./ (3 - rho).^3 .* (27 + 1.5 * (3 - rho - 3 ./ rho) .* (u.^2 + v.^2)  + E - 0.5 / tc.We * (w1.^2 + w2.^2)) - 6;
            case 2
                dp2 = 3 / tc.Cv * u .* (3 - 2 * rho) ./ rho ./ (3 - rho).^2;
            case 3
                dp2 = 3 / tc.Cv * v .* (3 - 2 * rho) ./ rho ./ (3 - rho).^2;
            case 4
                dp2 = 3 / tc.Cv ./ (3 - rho).^2;
            case 5
                dp2 = -3 / tc.Cv / tc.We * w1 ./ (3 - rho).^2;
            case 6
                dp2 = -3 / tc.Cv / tc.We * w2 ./ (3 - rho).^2;
            otherwise
                error('Wrong index')
        end
    case 2
        switch j
            case 1
                dp2 = 3 / tc.Cv * u .* (3 - 2 * rho) ./ rho ./ (3 - rho).^2;
            case 2
                dp2 = 3 / tc.Cv ./ rho ./ (rho - 3);
            case 3
                dp2 = zeros(size(rho));
            case 4
                dp2 = zeros(size(rho));
            case 5
                dp2 = zeros(size(rho));
            case 6
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end      
    case 3
        switch j
            case 1
                dp2 = 3 / tc.Cv * v .* (3 - 2 * rho) ./ rho ./ (3 - rho).^2;
            case 2
                dp2 = zeros(size(rho));
            case 3
                dp2 = 3 / tc.Cv ./ rho ./ (rho - 3);    
            case 4
                dp2 = zeros(size(rho));
            case 5
                dp2 = zeros(size(rho));
            case 6
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end
    case 4
        switch j
            case 1
                dp2 = 3 / tc.Cv ./ (3 - rho).^2;
            case 2
                dp2 = zeros(size(rho));
            case 3
                dp2 = zeros(size(rho));
            case 4
                dp2 = zeros(size(rho));
            case 5
                dp2 = zeros(size(rho));
            case 6
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end     
    case 5
        switch j
            case 1
                dp2 = -3 / tc.Cv / tc.We * w1 ./ (3 - rho).^2;
            case 2
                dp2 = zeros(size(rho));
            case 3
                dp2 = zeros(size(rho));
            case 4
                dp2 = zeros(size(rho));
            case 5
                dp2 = 3 / tc.Cv / tc.We ./ (rho - 3);
            case 6
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end  
    case 6
        switch j
            case 1
                dp2 = -3 / tc.Cv / tc.We * w2 ./ (3 - rho).^2;
            case 2
                dp2 = zeros(size(rho));
            case 3
                dp2 = zeros(size(rho));
            case 4
                dp2 = zeros(size(rho));
            case 5
                dp2 = zeros(size(rho));
            case 6
                dp2 = 3 / tc.Cv / tc.We ./ (rho - 3);
            otherwise
                error('Wrong index')
        end  
    otherwise
        error('Wrong index')
end

end