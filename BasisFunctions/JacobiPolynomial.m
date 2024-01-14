% Copyright (C) Xiangyi Meng.

function r = JacobiPolynomial(x, alpha, beta, poly_deg, der_order)

if(poly_deg < 0)
    error('Polynomial order of Jacobi polynomial must be >= 0')
end

if(poly_deg > 0)
    if(floor(poly_deg)~= poly_deg)
        error('Polynomial order of Jacobi polynomial must be an integer')
    end
end

if(alpha <= -1)
    error('First index of Jacobi polynomial must be larger than -1')
end

if(beta <= -1)
    error('Second index of Jacobi polynomial must be larger than -1')
end

switch der_order
    case 0
        if (poly_deg == 0)
            r = ones(size(x));
        elseif (poly_deg == 1)
            r = 0.5 * (alpha - beta + (alpha + beta + 2) * x);
        else
            n  = poly_deg - 1;
            An = 0.5 * (2 * n + alpha + beta + 1) * (2 * n + alpha + beta + 2) / ((n + 1) * (n + alpha + beta + 1));
            Bn = 0.5 * (alpha^2 - beta^2) * (2 * n + alpha + beta + 1) / ((n + 1) * (n + alpha + beta + 1) * (2 * n + alpha + beta));
            Cn = (n + alpha) * (n + beta) * (2 * n + alpha + beta + 2) / ((n + 1) * (n + alpha + beta + 1) * (2 * n + alpha + beta));
            
            r = (An * x + Bn) .* JacobiPolynomial(x, alpha, beta, n, der_order) - Cn * JacobiPolynomial(x, alpha, beta, n - 1, der_order);
        end
    case 1
        if (poly_deg == 0)
            r = zeros(size(x));
        else
            r = 0.5 * (poly_deg + alpha + beta + 1) * JacobiPolynomial(x, alpha + 1, beta + 1, poly_deg - 1, 0);
        end     
    case 2
        if (poly_deg == 0 || poly_deg == 1)
            r = zeros(size(x));
        else
            r = 0.25 * (poly_deg + alpha + beta + 1) * (poly_deg + alpha + beta + 2) * JacobiPolynomial(x, alpha + 2, beta + 2, poly_deg - 2, 0);
        end
    otherwise
        error('Not implemented derivative order of Jacobi polynomial')
end

end