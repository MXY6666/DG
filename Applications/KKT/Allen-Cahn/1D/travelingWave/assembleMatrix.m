function [IME, M, A1, S2, S] = assembleMatrix(msh, md, flux, tc, bs, mtol)

% mass matrix at reference line
ME = computeElementMatrix_refLine(0, 0, bs);

% inverse of mass matrix at reference line
IME = mychop(inv(ME), mtol);

% mass matrix
M = assembleElementMatrix_line(msh, 0, 0, bs, false, ME, mtol);

% inverse of mass matrix 
IM = assembleElementMatrix_line(msh, 0, 0, bs, true, IME, mtol);

% element contributions for stiffness matrix
SE = assembleElementMatrix_line(msh, 0, 1, bs, false, [], mtol);

% internal face contributions for stiffness matrix
SIFL = assembleIntFaceMatrix_line(msh, 'l', bs, md, [], mtol);
SIFR = assembleIntFaceMatrix_line(msh, 'r', bs, md, [], mtol);

SIF1 = (0.5 - flux.alpha) * SIFL + (0.5 + flux.alpha) * SIFR;
SIF2 = (0.5 + flux.alpha) * SIFL + (0.5 - flux.alpha) * SIFR;

% boundary face contributions for stiffness matrix
SBF = assembleBndFaceMatrix_line(msh, 'l', bs, md, [], 2, mtol);

% global matrices
A1 = mychop(IM * (SIF1 - SE), mtol);
S2 = tc.nu * mychop(SE - SIF2 - SBF, mtol);
S  = S2 * A1;

end