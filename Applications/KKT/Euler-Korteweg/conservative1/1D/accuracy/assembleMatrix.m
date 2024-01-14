function [A11, A12, Acap] = assembleMatrix(msh, md, flux, tc, bs, IME, mtol)

% inverse of mass matrix
IM = assembleElementMatrix_line(msh, 0, 0, bs, true, IME, mtol);

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% face contributions for stiffness matrix
SFL = assembleIntFaceMatrix_line(msh, 'l', bs, md, [], mtol);
SFR = assembleIntFaceMatrix_line(msh, 'r', bs, md, [], mtol);

SF11 = (0.5 - flux.alpha1) * SFL + (0.5 + flux.alpha1) * SFR;
SF12 = (0.5 + flux.alpha1) * SFL + (0.5 - flux.alpha1) * SFR;

SF21 = (0.5 - flux.alpha2) * SFL + (0.5 + flux.alpha2) * SFR;
SF22 = (0.5 + flux.alpha2) * SFL + (0.5 - flux.alpha2) * SFR;

% the global matrices
A11  = mychop(IM * (-SE + SF11), mtol);
A12  = mychop(IM * (-SE + SF12), mtol);
A21  = mychop(IM * (-SE + SF21), mtol);
A22  = mychop(IM * (-SE + SF22), mtol);
Acap = tc.kappa * mychop(A11 * A22 * A21, mtol);

end