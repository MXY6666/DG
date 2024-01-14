% Copyright (C) Xiangyi Meng.

function r = BernsteinPolynomial(x, n, k, der_order)

if(min(x) < 0 || max(x) > 1)
    error('Bernstein polynomials can only be computed for 0 <= x <= 1')
end

switch der_order
    case 0
        if(k < 0 || k > n)
            r = zeros(size(x));
            return
        end

        if(n == 0 && k == 0)
            r = ones(size(x));
        else
            r = x .* BernsteinPolynomial(x, n - 1, k - 1, 0) + (1 - x) .* BernsteinPolynomial(x, n - 1, k, 0);
        end        
    case 1
        r = n * (BernsteinPolynomial(x, n - 1, k - 1, 0) - BernsteinPolynomial(x, n - 1, k, 0));
    otherwise
        error('Not implemented derivative order of Jacobi polynomial')
end

end