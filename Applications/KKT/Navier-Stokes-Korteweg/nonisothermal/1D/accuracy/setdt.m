function dt = setdt(msh, md, U, W, t, CFL, tp, tc, bs)

% Compute reference wave velocity
waveVel = computeRefWaveVelocity(msh, md, U, W, tc, bs);

% Set the time step
dt = CFL * md.meshSize / waveVel;
if (t + dt) > tp
    dt = tp - t;
end

end

