function dt = setdt(msh, t, cfl, tp)

% Set the time step
dt = cfl * min(msh.elemLength(:, msh.LElems));

if (t + dt) > tp
    dt = tp - t;
end

end