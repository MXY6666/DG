function dt = setdt(msh, t, tc, cfl, tp)

% Set the time step
dt = cfl * min(msh.elemLength(:, msh.LElems)) / abs(tc.beta);

if (t + dt) > tp
    dt = tp - t;
end

end