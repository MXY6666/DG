function dt = setdt(msh, t, tc, cfl, tp)

% Set the time step
hmin = min(msh.elemLength(:, msh.LElems), [], 'all');
dt = cfl * hmin / max(abs(tc.beta));

if (t + dt) > tp
    dt = tp - t;
end

end