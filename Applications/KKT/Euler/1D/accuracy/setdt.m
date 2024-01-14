function dt = setdt(msh, md, U, t, CFL, tp, tc, bs)

% Compute reference wave velocity
waveVel = computeRefWaveVelocity(msh, md, U, tc, bs);

% Set the time step
dt = CFL * md.meshSize / waveVel;
if (t + dt) > tp
    dt = tp - t;
end

end

