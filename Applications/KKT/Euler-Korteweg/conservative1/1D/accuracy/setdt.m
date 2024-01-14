% Set the time step
function dt = setdt(msh, U, t, tc, quad1, bs, cfl, tp, isEx)

if isEx
    dt = cfl * min(msh.elemLength(:, msh.LElems)).^2;
else
    um = 0.5 * quad1.weights' * (bs.phi{1} * reshape(U, [bs.nb, 2 * msh.nLElems]));
    alpha = max(getEigenMax(um(1 : msh.nLElems), um(msh.nLElems + 1 : end), tc));
    dt = cfl * min(msh.elemLength(:, msh.LElems)) / ((alpha - 1) * (alpha > 1.0e-9) + 1);
end

if (t + dt) > tp
    dt = tp - t;
end

end