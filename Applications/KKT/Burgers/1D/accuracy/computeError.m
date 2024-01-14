function [L2err, Lierr] = computeError(msh, U, t, tc, quad, bs)

% Compute the exact solution at Gauss points in each element
ct = msh.elemCenter(:, msh.LElems);
h  = msh.elemLength(:, msh.LElems);
exact = zeros(quad.np, msh.nLElems);
for ie = 1 : msh.nLElems
    for j = 1 : quad.np
        exact(j, ie) = tc.u(ct(ie) + 0.5 * h(ie) * quad.points(j), t);
    end
end

% Compute the numerical solution at Gauss points in each element
numerical = bs.phi{1} * U;

% Compute L2 error
L2err = sqrt(sum(msh.elemJac(:, msh.LElems) .* (quad.weights' * (exact - numerical).^2)));

% Compute L-infinity error
Lierr = max(max(abs(exact - numerical))); 

end