function [AE, A2, A3, A32] = assembleMatrix(msh, md, flux, tc, bs, IME, mtol)

% inverse of mass matrix
IM = assembleElementMatrix_line(msh, 0, 0, bs, true, IME, mtol);

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% face contributions for stiffness matrix
SFL = assembleIntFaceMatrix_line(msh, 'l', bs, md, [], mtol);
SFR = assembleIntFaceMatrix_line(msh, 'r', bs, md, [], mtol);

SF21 = (0.5 - flux.alpha2) * SFL + (0.5 + flux.alpha2) * SFR;
SF22 = (0.5 + flux.alpha2) * SFL + (0.5 - flux.alpha2) * SFR;
SF31 = (0.5 - flux.alpha3) * SFL + (0.5 + flux.alpha3) * SFR;
SF32 = (0.5 + flux.alpha3) * SFL + (0.5 - flux.alpha3) * SFR;

% the global matrices
AE  = mychop(IM * (-SE), mtol);
A21 = mychop(IM * (-SE + SF21), mtol);
A22 = mychop(IM * (-SE + SF22), mtol);
A31 = mychop(IM * (-SE + SF31), mtol);
A32 = mychop(IM * (-SE + SF32), mtol);
A32 = 4 / 3 / tc.Re * A32;
A2  = 1 / tc.We * mychop(A22 * A21, mtol);
A3  = mychop(A32 * A31, mtol);

end