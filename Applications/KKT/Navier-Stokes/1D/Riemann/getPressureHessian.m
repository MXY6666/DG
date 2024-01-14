% Compute the Hessian matrix of pressure with respect to the conserved
% variables
function dp2 = getPressureHessian(rho, u, i, j, tc)

switch i
    case 1
        switch j
            case 1
                dp2 = (1 - tc.gamma) * u.^2 ./ rho;
            case 2
                dp2 = (tc.gamma - 1) * u ./ rho;
            case 3
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end
    case 2
        switch j
            case 1
                dp2 = (tc.gamma - 1) * u ./ rho;
            case 2
                dp2 = (1 - tc.gamma) ./ rho;
            case 3
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end        
    case 3
        switch j
            case 1
                dp2 = zeros(size(rho));
            case 2
                dp2 = zeros(size(rho));
            case 3
                dp2 = zeros(size(rho));
            otherwise
                error('Wrong index')
        end        
    otherwise
        error('Wrong index')
end

end