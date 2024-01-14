function U = getStageSolution(U0, KStage, dt, stage, bt)

if (stage.s == 0)
    U = U0;
elseif (bt.form == 1 && stage.s > 0)
    U = KStage(:, 1);
elseif (bt.form == 2 && stage.s > 0)
    U = U0 + dt * (KStage(:, 1 : stage.s) * bt.A(stage.s, 1 : stage.s)');
end

end