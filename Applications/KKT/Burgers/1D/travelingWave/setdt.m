function dt = setdt(msh, U, t, bs, cfl, tp)

% Determine the maximum absolute value of f'(uh) = uh on both sides of all faces
alpha = getEigenMax(msh, U, bs);

% Set the time step
dt = cfl * min(msh.elemLength(:, msh.LElems)) / ((alpha - 1) * (alpha > 1.0e-9) + 1);

if (t + dt) > tp
    dt = tp - t;
end

end